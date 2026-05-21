"""
main.py
Authors:
    Filipe Galiforni Silva (TU Delft)
    Ad Reniers (TU Delft)
    Sierd de Vries (TU Delft)
    Floris van Rees (Deltares)
    Caroline Hallin (Lund University).

TU Delft, 2025

Description:
    1D Version of the CS Model (Hallin et. al., 2019) adapted to include vegetation and hard structures
    to simulate hybrid dune-dike systems. The model simulates cross-shore profile evolution under the influence
    of waves, water levels, vegetation growth, and the presence of hard structures.

Usage (more info on the README file):
    python main.py path/to/config.yaml

Inputs (within the YAML configuration file):
    - Initial cross-shore profile data
    - Wave and water level data
    - Vegetation parameters
    - Hard structure parameters

Outputs:
    - Updated cross-shore profile over time
    - Updated vegetation coverage over time

Dependencies:
    - numpy
    - scipy
    - matplotlib
    - shapely
    - pandas
    - netCDF4

License:
    Copyright (c) 2025, Delft University of Technology.
    All rights reserved.

References:
HALLIN, C., LARSON, M. & HANSON, H. 2019. Simulating beach and dune evolution at decadal to centennial scale under rising sea levels. PLoS One, 14, e0215651. doi: 10.1371/journal.pone.0215651

"""

############################################################################
# IMPORTS
############################################################################

import sys
import os
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.colors import ListedColormap
from shapely.geometry import LineString
from scipy.interpolate import interp1d
import pandas as pd
import netCDF4 as nc
import yaml
import zlib
from types import SimpleNamespace


############################################################################
# PATH SETUP
############################################################################

# if len(sys.argv) < 2:
#     raise RuntimeError(
#         "Usage: python main.py path/to/config.yaml"
#     )

#config_path = sys.argv[1]

package_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(package_dir)  # parent folder

if project_root not in sys.path:
    sys.path.insert(0, project_root)

from cs_model.functions import *
from cs_model.mass_conservation import calculate_compartment_volumes, calculate_timestep_conservation
import pickle

# Persistent vegetation rendering state used by snapshot plotting so stems
# evolve smoothly from one timestep to the next within the same run.
_VEG_RENDER_STATE = {}

############################################################################
# INITIALIZATION OF PARAMETERS
############################################################################

def load_config(config_path):
    """
    Load YAML configuration file from a given path.
    """
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    return cfg


def _save_profile_snapshot(dune_x_init, dune_z_init, dune_x_curr, dune_z_curr,
                           substrate_x, substrate_z, output_base, timestep,
                           progress_pct, bc_info=None, vegetation_info=None,
                           y_lim_top=None, y_lim_bottom=None,
                           x_lim_right=None, x_lim_left=None):
    """
    Save a snapshot plot of the current cross-shore profile vs the initial profile.

    Parameters
    ----------
    dune_x_init : dict
        Cross-shore positions of dune pivots at the start of the simulation.
    dune_z_init : dict
        Elevations of dune pivots at the start of the simulation.
    dune_x_curr : dict
        Cross-shore positions of dune pivots at the current timestep.
    dune_z_curr : dict
        Elevations of dune pivots at the current timestep.
    substrate_x : list
        Cross-shore positions of the non-erodible substrate.
    substrate_z : list
        Elevations of the non-erodible substrate.
    output_base : str
        Base path for output files (without extension).
    timestep : int
        Completed timestep count (1-based; used directly in the filename so
        filenames are human-readable, e.g. ``output_5000.png``).
    progress_pct : float
        Current simulation progress in percent (used in labels).
    bc_info : dict, optional
        Boundary condition values to annotate on the plot.  Expected keys
        (all optional): ``wave_height`` (Hs, m), ``wave_period`` (Tp, s),
        ``runup`` (m), ``total_water_level`` (TWL, m).
        When ``total_water_level`` is present a water-level fill is drawn
        on the seaward side of the profile using the sea colour.
    vegetation_info : dict, optional
        Normalized vegetation coverage values (0–1) for three dune zones:
        ``rhov_df`` (dune front/ramp), ``rhov_dt`` (dune top), and
        ``rhov_db`` (dune back/lee).  Optional second-species values
        (``rhov_df2``, ``rhov_dt2``, ``rhov_db2``) are summed and clipped to 1.
        Rendered as short curved stems anchored to the current dune profile.
    y_lim_top : float or None, optional
        Fixed upper limit for the y-axis.  When provided, all snapshots share
        the same vertical range so that sudden water-level events do not cause
        the axis to rescale between frames.  When ``None`` (default), the
        upper limit is determined automatically by matplotlib.
    y_lim_bottom : float or None, optional
        Fixed lower limit for the y-axis.  When provided, every snapshot uses
        this value.  When ``None`` (default), 0 is used.
    x_lim_right : float or None, optional
        Fixed right limit for the x-axis.  When provided, every snapshot uses
        this value so GIF frames have consistent horizontal extent.  When
        ``None`` (default) the limit is computed dynamically from the current
        shoreline position.
    x_lim_left : float or None, optional
        Fixed left limit for the x-axis.  When provided, every snapshot uses
        this value.  When ``None`` (default), 0 is used.

    Returns
    -------
    plot_filename : str or None
        Path of the saved PNG file, or ``None`` if the snapshot was skipped.
    """
    c_init  = "#1f77b4"   # blue
    c_final = "#d62728"   # red
    c_sub   = "0.3"       # grey
    c_fill_pos = "#ff7f0e"  # orange (deposition)
    c_fill_neg = "#7b2cbf"  # purple (erosion)
    c_sea   = "#b7dff5"   # soft ocean blue
    c_sand  = "#f4e6c8"   # very light sand
    
    _profile_order = [
        'back_ref', 'lee_foot', 'lee_crest', 'crest_temp',
        'seaward_crest', 'ramp_dike_temp', 'ramp_crest',
        'ramp_toe_temp', 'ramp_toe', 'shoreline'
    ]

    def _ordered_profile_points(dune_x, dune_z, start_key=None):
        keys = [
            key for key in _profile_order
            if key in dune_x and key in dune_z
            and dune_x[key] is not None and dune_z[key] is not None
        ]
        if start_key in keys:
            keys = keys[keys.index(start_key):]
        x = np.array([dune_x[key] for key in keys], dtype=float)
        z = np.array([dune_z[key] for key in keys], dtype=float)
        if x.size == 0:
            return x, z

        # Drop any coordinate that is NaN or Inf – a NaN breaks the plotted line.
        _finite = np.isfinite(x) & np.isfinite(z)
        keys  = [k for k, ok in zip(keys, _finite) if ok]
        x, z  = x[_finite], z[_finite]
        if x.size == 0:
            return x, z

        _present = {key: idx for idx, key in enumerate(keys)}
        if 'crest_temp' in _present and 'seaward_crest' in _present:
            ict, isc = _present['crest_temp'], _present['seaward_crest']
            if x[ict] > x[isc]:
                x[ict] = x[isc]
        if 'seaward_crest' in _present and 'ramp_dike_temp' in _present:
            isc, irdt = _present['seaward_crest'], _present['ramp_dike_temp']
            if x[isc] > x[irdt]:
                x[isc] = x[irdt]
        if 'ramp_dike_temp' in _present and 'ramp_crest' in _present:
            irdt, irc = _present['ramp_dike_temp'], _present['ramp_crest']
            if x[irc] < x[irdt]:
                x[irc] = x[irdt]

        order = np.argsort(x, kind='stable')
        return x[order], z[order]

    def _fill_profile_points(dune_x, dune_z, start_key=None):
        x, z = _ordered_profile_points(dune_x, dune_z, start_key=start_key)
        if x.size == 0:
            return x, z, x, z
        if x.size == 1:
            return x, z, x, z

        n_per_segment = 100
        x_dense = np.concatenate([
            np.linspace(x[i], x[i + 1], n_per_segment, endpoint=False)
            for i in range(len(x) - 1)
        ] + [[x[-1]]])
        z_dense = np.concatenate([
            np.linspace(z[i], z[i + 1], n_per_segment, endpoint=False)
            for i in range(len(z) - 1)
        ] + [[z[-1]]])

        x_unique, inverse = np.unique(x_dense, return_inverse=True)
        z_upper = np.full(x_unique.shape, -np.inf, dtype=float)
        np.maximum.at(z_upper, inverse, z_dense)
        return x, z, x_unique, z_upper

    x_init, z_init, x_init_fill, z_init_fill = _fill_profile_points(dune_x_init, dune_z_init)
    x_curr, z_curr, x_curr_fill, z_curr_fill = _fill_profile_points(dune_x_curr, dune_z_curr)
    _, _, x_init_body_fill, z_init_body_fill = _fill_profile_points(
        dune_x_init, dune_z_init, start_key='lee_foot'
    )
    _, _, x_curr_body_fill, z_curr_body_fill = _fill_profile_points(
        dune_x_curr, dune_z_curr, start_key='lee_foot'
    )

    if x_init_fill.size == 0 or x_curr_fill.size == 0:
        print(f"Profile snapshot at {progress_pct:.0f}% skipped (profiles are empty or degenerate).")
        return None

    xg_start = min(x_init_fill.min(), x_curr_fill.min())
    xg_end   = max(x_init_fill.max(), x_curr_fill.max())
    if x_init_body_fill.size == 0 or x_curr_body_fill.size == 0 or xg_start >= xg_end:
        # Profiles are empty or degenerate – skip this snapshot
        print(f"Profile snapshot at {progress_pct:.0f}% skipped (profiles are empty or degenerate).")
        return None

    xg = np.linspace(xg_start, xg_end, 500)

    #z_init_g = np.interp(xg, x_init, z_init)
    #z_curr_g = np.interp(xg, x_curr, z_curr)
    z_sub_g = np.interp(xg, substrate_x, substrate_z)

    # Full-profile fills should reach the latest visible profile everywhere,
    # including vertical/duplicate-x sections and the landward back-ref segment.
    z_init_fill_g = np.interp(xg, x_init_fill, z_init_fill, left=np.nan, right=np.nan)
    z_curr_fill_g = np.interp(xg, x_curr_fill, z_curr_fill, left=np.nan, right=np.nan)

    # Change fills stay bounded from lee_foot onward so the back-ref segment does
    # not create spurious erosion/deposition shading landward of the dune body.
    z_init_g = np.interp(xg, x_init_body_fill, z_init_body_fill, left=np.nan, right=np.nan)
    z_curr_g = np.interp(xg, x_curr_body_fill, z_curr_body_fill, left=np.nan, right=np.nan)
    
    # Masks for robust change filling
    init_valid = np.isfinite(z_init_g)
    curr_valid = np.isfinite(z_curr_g)
    both_valid = init_valid & curr_valid

    fig, ax = plt.subplots(figsize=(10, 4.8))
    
    # Light sand body under the final profile
    ax.fill_between(xg,
                0.0,
                z_curr_fill_g,
                where=np.isfinite(z_curr_fill_g),
                facecolor=c_sand,
                edgecolor='none',
                alpha=0.6,
                zorder=0.5)

    ax.fill_between(xg, z_init_g, z_curr_g,
                    #where=(z_curr_g > z_init_g),
                    where=(both_valid & (z_curr_g > z_init_g)),
                    color=c_fill_pos,
                    interpolate=True, alpha=0.35, zorder=1, label="Deposition")

    ax.fill_between(xg, z_init_g, z_curr_g,
                    #where=(z_curr_g < z_init_g),
                    where=(both_valid & (z_curr_g < z_init_g)),
                    color=c_fill_neg,
                    interpolate=True, alpha=0.35, zorder=1, label="Erosion")

    #ax.plot(x_init[1:], z_init[1:], color=c_init, lw=2.8, zorder=3, label="Initial profile")
    #ax.plot(x_curr[1:], z_curr[1:], color=c_final, lw=2.8, zorder=3,
    #        label=f"Profile at {progress_pct:.0f}%")
  # Deposition in extension zones where only final profile exists (progradation/growth)
    ax.fill_between(xg, z_sub_g, z_curr_g,
                    where=(curr_valid & ~init_valid & (z_curr_g > z_sub_g)),
                    color=c_fill_pos,
                    interpolate=True, alpha=0.35, zorder=1)

    # Fallback deposition in final-only zones using edge level of initial profile.
    # This captures shoreline progradation beyond initial x-limits when substrate
    # capping would otherwise suppress the fill.
    z_init_edge = np.where(xg < x_init_body_fill.min(), z_init_body_fill[0], z_init_body_fill[-1])
    ax.fill_between(xg, z_init_edge, z_curr_g,
                    where=(curr_valid & ~init_valid & (z_curr_g > z_init_edge) & ~(z_curr_g > z_sub_g)),
                    color=c_fill_pos,
                    interpolate=True, alpha=0.35, zorder=1)

    # Erosion in extension zones where only initial profile exists (retreat)
    ax.fill_between(xg, z_sub_g, z_init_g,
                    where=(init_valid & ~curr_valid & (z_init_g > z_sub_g)),
                    color=c_fill_neg,
                    interpolate=True, alpha=0.35, zorder=1)

    ax.plot(x_init, z_init, color=c_init, lw=2.8, zorder=3, label="Initial profile")
    ax.plot(x_curr, z_curr, color=c_final, lw=2.8, zorder=3,
            label=f"Profile at {progress_pct:.0f}%")
    
    # ---- Vegetation overlay ----
    if vegetation_info is None:
        vegetation_info = {}

    _veg_cov = {
        'df': float(
            np.clip(
                vegetation_info.get('rhov_df', 0.0) + vegetation_info.get('rhov_df2', 0.0),
                0.0, 1.0
            )
        ),
        'dt': float(
            np.clip(
                vegetation_info.get('rhov_dt', 0.0) + vegetation_info.get('rhov_dt2', 0.0),
                0.0, 1.0
            )
        ),
        'db': float(
            np.clip(
                vegetation_info.get('rhov_db', 0.0) + vegetation_info.get('rhov_db2', 0.0),
                0.0, 1.0
            )
        ),
    }

    _veg_cmap = ListedColormap(plt.cm.Greens(np.linspace(0.22, 0.98, 256)))
    _veg_norm = plt.Normalize(vmin=0.0, vmax=1.0)
    _any_veg_plotted = False
    _veg_run_key = os.path.abspath(str(output_base))
    _veg_run_state = _VEG_RENDER_STATE.setdefault(_veg_run_key, {})

    def _plot_veg_band(x0, x1, cov, band_key):
        nonlocal _any_veg_plotted
        _COV_ZERO_EPS = 1e-6
        _COV_PLOT_MIN = 1e-3

        def _clear_band_state(_state):
            _state["eff_cov"] = 0.0
            _state["stems_u"] = np.array([], dtype=float)
            _state["bend_lat"] = np.array([], dtype=float)
            _state["tip_lat"] = np.array([], dtype=float)
            _state["h_lat"] = np.array([], dtype=float)
            _state["cluster_centers"] = np.array([], dtype=float)
            _state["cluster_ids"] = np.array([], dtype=int)

        _band_state = _veg_run_state.get(band_key, {})
        if not np.isfinite(cov) or cov <= _COV_ZERO_EPS:
            _clear_band_state(_band_state)
            _veg_run_state[band_key] = _band_state
            return
        xa, xb = float(min(x0, x1)), float(max(x0, x1))
        if xb - xa <= 1e-9:
            return

        _seed = int(zlib.crc32(f"{_veg_run_key}:{band_key}".encode("utf-8")) & 0xFFFFFFFF)
        _veg_rng = _band_state.get("rng")
        if _veg_rng is None:
            _veg_rng = np.random.default_rng(_seed)
            _band_state["rng"] = _veg_rng

        _eff_cov_prev = float(_band_state.get("eff_cov", cov))
        _eff_cov = float(np.clip(0.82 * _eff_cov_prev + 0.18 * cov, 0.0, 1.0))
        if _eff_cov < _COV_PLOT_MIN:
            _clear_band_state(_band_state)
            _veg_run_state[band_key] = _band_state
            return
        _band_state["eff_cov"] = _eff_cov

        _target_stems = int(
            np.clip(
                np.round((xb - xa) * (2.20 * (_eff_cov ** 1.30))),
                0, 160
            )
        )
        if _target_stems <= 0:
            _clear_band_state(_band_state)
            _veg_run_state[band_key] = _band_state
            return
        _MAX_PATCH = 6  # maximum stems per vegetation patch
        _PATCH_SIG = 0.013 + 0.007 * _eff_cov  # spatial spread of stems within a patch

        stems_u = np.asarray(_band_state.get("stems_u", np.array([])), dtype=float)
        bend_lat = np.asarray(_band_state.get("bend_lat", np.array([])), dtype=float)
        tip_lat = np.asarray(_band_state.get("tip_lat", np.array([])), dtype=float)
        h_lat = np.asarray(_band_state.get("h_lat", np.array([])), dtype=float)
        cluster_centers = np.asarray(_band_state.get("cluster_centers", np.array([])), dtype=float)
        cluster_ids = np.asarray(_band_state.get("cluster_ids", np.array([])), dtype=int)

        if stems_u.size == 0:
            # Need at least enough clusters to hold all stems at max _MAX_PATCH each.
            _min_clusters = int(np.ceil(_target_stems / _MAX_PATCH))
            _n_cluster = max(
                _min_clusters,
                int(np.clip(np.round(2 + (xb - xa) * (0.06 + 0.14 * _eff_cov)), _min_clusters, 24)),
            )
            # Evenly spaced cluster centers with small jitter.
            _base = np.linspace(0.08, 0.92, _n_cluster)
            _jitter = _veg_rng.uniform(-0.04, 0.04, size=_n_cluster)
            cluster_centers = np.clip(_base + _jitter, 0.02, 0.98)
            # Distribute stems across clusters: round-robin so each gets ≤ _MAX_PATCH.
            _cluster_counts = np.zeros(_n_cluster, dtype=int)
            for _si in range(_target_stems):
                _ci = _si % _n_cluster
                _cluster_counts[_ci] += 1
            _all_u, _all_cid = [], []
            for _ci, (_ctr, _cnt) in enumerate(zip(cluster_centers, _cluster_counts)):
                if _cnt <= 0:
                    continue
                _u = np.clip(_veg_rng.normal(_ctr, _PATCH_SIG, size=_cnt), 0.0, 1.0)
                _all_u.append(_u)
                _all_cid.extend([_ci] * _cnt)
            stems_u = np.concatenate(_all_u) if _all_u else cluster_centers[:1].copy()
            cluster_ids = np.array(_all_cid, dtype=int)
            bend_lat = _veg_rng.normal(0.0, 1.0, size=stems_u.size)
            tip_lat = _veg_rng.normal(0.0, 1.0, size=stems_u.size)
            h_lat = _veg_rng.normal(0.0, 1.0, size=stems_u.size)
        else:
            # Drift stem positions and shape attributes.
            stems_u = np.clip(stems_u + _veg_rng.normal(0.0, 0.008, size=stems_u.size), 0.0, 1.0)
            bend_lat = 0.90 * bend_lat + 0.10 * _veg_rng.normal(0.0, 1.0, size=bend_lat.size)
            tip_lat = 0.90 * tip_lat + 0.10 * _veg_rng.normal(0.0, 1.0, size=tip_lat.size)
            h_lat = 0.88 * h_lat + 0.12 * _veg_rng.normal(0.0, 1.0, size=h_lat.size)

            # Rebuild cluster_ids if missing (backwards compat or first update).
            if cluster_centers.size == 0:
                cluster_centers = stems_u.copy()
                cluster_ids = np.arange(stems_u.size, dtype=int)
            _n_cluster = cluster_centers.size

            if stems_u.size < _target_stems:
                _n_add = _target_stems - stems_u.size
                _counts_per = np.bincount(cluster_ids, minlength=_n_cluster)
                _capacity = np.maximum(0, _MAX_PATCH - _counts_per)
                # If existing clusters are all full, spawn new ones.
                if _capacity.sum() < _n_add:
                    _extra = int(np.ceil((_n_add - _capacity.sum()) / _MAX_PATCH))
                    _new_ctrs = _veg_rng.uniform(0.05, 0.95, size=_extra)
                    cluster_centers = np.concatenate([cluster_centers, _new_ctrs])
                    _capacity = np.concatenate([_capacity, np.full(_extra, _MAX_PATCH, dtype=int)])
                    _counts_per = np.concatenate([_counts_per, np.zeros(_extra, dtype=int)])
                    _n_cluster += _extra
                # Add to least-filled clusters first.
                _add_counts = np.zeros(_n_cluster, dtype=int)
                _remaining = _n_add
                for _ci in np.argsort(_counts_per):
                    if _remaining <= 0:
                        break
                    _give = min(int(_capacity[_ci]), _remaining)
                    _add_counts[_ci] = _give
                    _remaining -= _give
                _new_u_list, _new_cid_list = [], []
                for _ci, _cnt in enumerate(_add_counts):
                    if _cnt <= 0:
                        continue
                    _u = np.clip(_veg_rng.normal(cluster_centers[_ci], _PATCH_SIG, size=_cnt), 0.0, 1.0)
                    _new_u_list.append(_u)
                    _new_cid_list.extend([_ci] * _cnt)
                if _new_u_list:
                    _new_u = np.concatenate(_new_u_list)
                    _new_cid = np.array(_new_cid_list, dtype=int)
                    stems_u = np.concatenate([stems_u, _new_u])
                    cluster_ids = np.concatenate([cluster_ids, _new_cid])
                    bend_lat = np.concatenate([bend_lat, _veg_rng.normal(0.0, 1.0, size=_new_u.size)])
                    tip_lat = np.concatenate([tip_lat, _veg_rng.normal(0.0, 1.0, size=_new_u.size)])
                    h_lat = np.concatenate([h_lat, _veg_rng.normal(0.0, 1.0, size=_new_u.size)])

            elif stems_u.size > _target_stems:
                _n_remove = stems_u.size - _target_stems
                _counts_per = np.bincount(cluster_ids, minlength=_n_cluster)
                # Remove from most-populated clusters first.
                _remove_counts = np.zeros(_n_cluster, dtype=int)
                _remaining = _n_remove
                for _ci in np.argsort(-_counts_per):
                    if _remaining <= 0:
                        break
                    _take = min(int(_counts_per[_ci]), _remaining)
                    _remove_counts[_ci] = _take
                    _remaining -= _take
                _keep_mask = np.ones(stems_u.size, dtype=bool)
                for _ci, _cnt in enumerate(_remove_counts):
                    if _cnt <= 0:
                        continue
                    _in_cluster = np.where(cluster_ids == _ci)[0]
                    _keep_mask[_in_cluster[:_cnt]] = False
                stems_u = stems_u[_keep_mask]
                cluster_ids = cluster_ids[_keep_mask]
                bend_lat = bend_lat[_keep_mask]
                tip_lat = tip_lat[_keep_mask]
                h_lat = h_lat[_keep_mask]

        veg_color = _veg_cmap(_veg_norm(cov))
        veg_lw = 0.42 + 1.55 * cov
        veg_h = 0.33 + 1.65 * cov
        _bend_scale = max(0.70, 0.110 * (xb - xa)) * (0.90 + 1.30 * _eff_cov)
        _tip_scale = max(0.15, 0.035 * (xb - xa)) * (0.70 + 0.95 * _eff_cov)
        _plotted_in_band = False

        xv_stem = xa + stems_u * (xb - xa)
        zv_stem = np.interp(xv_stem, x_curr_fill, z_curr_fill, left=np.nan, right=np.nan)
        z_sub_stem = np.interp(xv_stem, substrate_x, substrate_z)
        _ok_stem = np.isfinite(zv_stem) & (zv_stem >= (z_sub_stem - 3e-2))

        for _x, _z, _b, _tj, _hl, _ok in zip(xv_stem, zv_stem, bend_lat, tip_lat, h_lat, _ok_stem):
            if not _ok:
                continue
            _h = veg_h * np.clip(1.0 + 0.22 * _hl, 0.52, 1.62)
            _curv = float(np.clip(_b, -2.0, 2.0)) * _bend_scale
            _tip_j = float(np.clip(_tj, -2.0, 2.0)) * _tip_scale
            _p0x, _p0z = _x, _z
            _p1x, _p1z = _x + _curv, _z + 0.56 * _h
            _p2x, _p2z = _x + 0.52 * _curv + _tip_j, _z + _h
            _t = np.linspace(0.0, 1.0, 9)
            _omt = 1.0 - _t
            _xc = (_omt * _omt) * _p0x + 2.0 * _omt * _t * _p1x + (_t * _t) * _p2x
            _zc = (_omt * _omt) * _p0z + 2.0 * _omt * _t * _p1z + (_t * _t) * _p2z
            ax.plot(_xc, _zc, color=veg_color, linewidth=veg_lw, alpha=0.95, zorder=6.2)
            _plotted_in_band = True

        # fallback to ensure visible stems in sparse/edge cases
        if not _plotted_in_band:
            _n_force = int(np.clip(np.round(2 + 6 * _eff_cov), 2, 10))
            _x_force = np.linspace(xa, xb, _n_force)
            _z_force = np.interp(_x_force, x_curr_fill, z_curr_fill, left=np.nan, right=np.nan)
            for _x, _z in zip(_x_force, _z_force):
                if not np.isfinite(_z):
                    continue
                _h = veg_h * np.clip(1.0 + _veg_rng.normal(0.0, 0.20), 0.60, 1.35)
                _curv = _veg_rng.uniform(-1.0, 1.0) * (1.20 * _bend_scale)
                _tip_j = _veg_rng.uniform(-1.0, 1.0) * (0.95 * _tip_scale)
                _p0x, _p0z = _x, _z
                _p1x, _p1z = _x + _curv, _z + 0.54 * _h
                _p2x, _p2z = _x + 0.5 * _curv + _tip_j, _z + _h
                _t = np.linspace(0.0, 1.0, 9)
                _omt = 1.0 - _t
                _xc = (_omt * _omt) * _p0x + 2.0 * _omt * _t * _p1x + (_t * _t) * _p2x
                _zc = (_omt * _omt) * _p0z + 2.0 * _omt * _t * _p1z + (_t * _t) * _p2z
                ax.plot(_xc, _zc, color=veg_color, linewidth=veg_lw, alpha=0.95, zorder=6.2)
            _plotted_in_band = True

        _band_state["stems_u"] = stems_u
        _band_state["bend_lat"] = bend_lat
        _band_state["tip_lat"] = tip_lat
        _band_state["h_lat"] = h_lat
        _band_state["cluster_centers"] = cluster_centers
        _band_state["cluster_ids"] = cluster_ids
        _veg_run_state[band_key] = _band_state
        _any_veg_plotted = _any_veg_plotted or _plotted_in_band

    def _first_finite(*vals):
        for _v in vals:
            if _v is not None and np.isfinite(_v):
                return float(_v)
        return np.nan

    _plot_veg_band(
        dune_x_curr.get('ramp_toe', np.nan),
        dune_x_curr.get('seaward_crest', dune_x_curr.get('ramp_crest', np.nan)),
        _veg_cov['df'],
        'df',
    )
    _plot_veg_band(
        dune_x_curr.get('seaward_crest', dune_x_curr.get('ramp_crest', np.nan)),
        dune_x_curr.get('lee_crest', np.nan),
        _veg_cov['dt'],
        'dt',
    )
    _plot_veg_band(
        _first_finite(dune_x_curr.get('lee_crest'), dune_x_curr.get('crest_temp')),
        _first_finite(dune_x_curr.get('lee_foot'), dune_x_curr.get('back_ref')),
        _veg_cov['db'],
        'db',
    )

    # ---- Boundary conditions ----
    # Use surface for water-level fills.
    # * Inside the current profile extent   → profile elevation.
    # * Seaward of the current shoreline    → sea level (0.0).  Using the
    #   substrate fallback here can produce elevated values in the ocean area
    #   (substrate extrapolation) which create a blank gap between the sea fill
    #   and the profile line.
    # * Landward NaN region (dike interior) → substrate elevation.
    _x_curr_shore = float(
        dune_x_curr.get(
            'shoreline',
            float(x_curr_fill[-1]) if x_curr_fill.size > 0 else xg_end,
        )
    )
    z_curr_for_sea = np.where(
        np.isfinite(z_curr_fill_g),
        z_curr_fill_g,
        np.where(xg > _x_curr_shore, 0.0, z_sub_g),
    )

    # Compute seaward extent of the sea fill (used for xlim regardless of bc_info)
    _x_ramp_toe_lim  = dune_x_curr.get('ramp_toe', x_curr_fill.max())
    _x_shoreline_lim = dune_x_curr.get('shoreline', x_curr_fill.max())
    x_sea_end = _x_shoreline_lim + 0.1 * (_x_shoreline_lim - _x_ramp_toe_lim)
    x_fill_right = x_lim_right if x_lim_right is not None else x_sea_end
    # Keep sea-fill x coordinates monotonic (required by fill_between) even if a
    # fixed shared x-limit is smaller than the current profile extent.
    x_fill_right = max(float(x_fill_right), float(xg[-1]) + 1e-9)

    if bc_info is not None:
        twl = bc_info.get('total_water_level')
        hs  = bc_info.get('wave_height')
        tp  = bc_info.get('wave_period')
        r   = bc_info.get('runup')
        wind_speed  = bc_info.get('wind_speed')
        wind_dir    = bc_info.get('wind_direction')

        x_ramp_toe  = dune_x_curr.get('ramp_toe', x_curr_fill.max())
        x_shoreline = dune_x_curr.get('shoreline', x_curr_fill.max())

        # Extended profile grid to include the seaward visual margin
        _xfill = np.concatenate([xg, [x_fill_right]])
        #_zfill = np.concatenate([z_curr_g, [0.0]])
        _zfill = np.concatenate([z_curr_for_sea, [0.0]])

        def _find_touch_x(xarr, zarr, level):
            """Return the seaward-most x where the profile first rises to *level*."""
            above = zarr >= level
            if not np.any(above):
                return xarr[0]
            last = int(np.where(above)[0][-1])
            if last < len(xarr) - 1:
                xa, xb = xarr[last], xarr[last + 1]
                za, zb = zarr[last], zarr[last + 1]
                if abs(zb - za) > 1e-10:
                    return xa + (level - za) / (zb - za) * (xb - xa)
            return xarr[last]

        _is_overwash = False
        if twl is not None:
            _dune_crest = dune_z_curr.get('seaward_crest')
            _is_overwash = (_dune_crest is not None and twl >= _dune_crest)
            if _is_overwash:
                # Overwash: fill from x=0 all the way seaward
                _z_back = dune_z_curr.get('back_ref',
                                          z_curr_for_sea[0] if len(z_curr_for_sea) > 0 else 0.0)
                _xfill_ow = np.concatenate([[0.0], xg, [x_fill_right]])
                _zfill_ow = np.concatenate([[_z_back], z_curr_for_sea, [0.0]])
                ax.fill_between(_xfill_ow, twl, _zfill_ow,
                                where=(_zfill_ow <= twl),
                                facecolor=c_sea, alpha=0.50, interpolate=True,
                                zorder=1, label=f"TWL = {twl:.2f} m (overwash)")
                ax.plot([0.0, x_fill_right], [twl, twl],
                        color=c_sea, lw=1.8, linestyle=':', zorder=5)
            else:
                # Find touch point then restrict fill to the seaward side only.
                _tx_twl = _find_touch_x(_xfill, _zfill, twl)
                # Insert the exact touch-point so fill_between starts precisely
                # at the profile/TWL intersection with no discretisation gap.
                _ins = int(np.searchsorted(_xfill, _tx_twl))
                _xfill_i = np.insert(_xfill, _ins, _tx_twl)
                _zfill_i = np.insert(_zfill, _ins, twl)
                ax.fill_between(_xfill_i, twl, _zfill_i,
                                where=(_xfill_i >= _tx_twl),
                                facecolor=c_sea, alpha=0.50, interpolate=True,
                                zorder=1, label=f"TWL = {twl:.2f} m")
                ax.plot([_tx_twl, x_fill_right], [twl, twl],
                        color=c_sea, lw=1.8, linestyle=':', zorder=5)

        # Overwash flow arrows (high-water events), aligned with seaward→landward flow.
        if twl is not None and twl > 5.0:
            _flow_y = max(0.2, twl - 0.35)
            if _is_overwash:
                _flow_x = np.array([0.78, 0.58, 0.38, 0.18]) * x_fill_right
                _flow_dx = np.full_like(_flow_x, -0.10 * x_fill_right)
            else:
                _flow_x = np.array([0.92, 0.80]) * x_fill_right
                _flow_dx = np.full_like(_flow_x, -0.08 * x_fill_right)
            ax.quiver(
                _flow_x, np.full_like(_flow_x, _flow_y),
                _flow_dx, np.zeros_like(_flow_dx),
                angles="xy", scale_units="xy", scale=1,
                width=0.0035, headwidth=5.2, headlength=6.8, headaxislength=5.8,
                color="#1f78b4", zorder=7,
            )
            _flow_label_x = 0.5 * (_flow_x[0] + _flow_x[-1]) + 0.5 * abs(_flow_dx[0])
            ax.text(
                _flow_label_x, _flow_y + 0.35, r"$\vec{U}_w$",
                ha="center", va="bottom", fontsize=11, color="#1f78b4", zorder=8,
            )

        # BC annotation box (upper-right)
        bc_lines = []
        if hs         is not None: bc_lines.append(f"Hs   = {hs:.2f} m")
        if tp         is not None: bc_lines.append(f"Tp   = {tp:.1f} s")
        if wind_speed is not None: bc_lines.append(f"U10  = {wind_speed:.1f} m/s")
        if wind_dir   is not None: bc_lines.append(f"Wdir = {wind_dir:.0f}\u00b0")
        if r          is not None: bc_lines.append(f"R    = {r:.2f} m")
        if twl        is not None: bc_lines.append(f"TWL  = {twl:.2f} m")

        if bc_lines:
            ax.text(0.98, 0.97, "\n".join(bc_lines),
                    transform=ax.transAxes,
                    fontsize=10.5,
                    verticalalignment='top',
                    horizontalalignment='right',
                    bbox=dict(boxstyle="round,pad=0.4",
                              facecolor="white",
                              edgecolor="0.4",
                              alpha=0.92))

        # Wind arrow: skip during high-water events to avoid overlap with flow arrows.
        if wind_speed is not None and (twl is None or twl <= 5.0):
            _xw0 = 500.0
            _yw0 = 6.0
            _wind_mag = 0.12 * x_fill_right
            ax.quiver(
                [_xw0], [_yw0], [-_wind_mag], [0.0],
                angles="xy", scale_units="xy", scale=1,
                width=0.0046, headwidth=6.2, headlength=7.8, headaxislength=6.8,
                color="0.2", zorder=9,
            )
            ax.text(
                _xw0 - 0.5 * _wind_mag, _yw0 + 0.65, r"$\vec{U}_{wind}$",
                ha="center", va="bottom", fontsize=11, color="0.2", zorder=10,
            )

    # Fill the substrate body with a solid colour to mask any erosion/deposition
    # fills or water fills that extend into the non-erodible dike structure.
    # Only fill up to the dike seaward toe (substrate_x[4] / substrate_z[4]);
    # the flat beach seaward of the dike should not be covered.
    # zorder=5 ensures the fill renders on top of the dune profile lines
    # (zorder=3) so that any dune geometry crossing inside the dike is hidden.
    _sub_z_min = min(z_init.min(), z_curr.min(), min(substrate_z)) - 1.0
    ax.fill_between(substrate_x[0:5], substrate_z[0:5], _sub_z_min,
                    facecolor='0.85', linewidth=0, zorder=5)
    # Redraw the substrate outline on top so it is always visible
    ax.plot(substrate_x[0:5], substrate_z[0:5],
            linestyle="--", color=c_sub, lw=2.0, zorder=6, label="Non-erodible substrate")

    ax.set_xlabel("Cross-shore distance (m)", fontsize=15)
    ax.set_ylabel("Elevation (m)", fontsize=15)
    ax.set_title(f"Cross-shore profile evolution – {progress_pct:.0f}% of simulation",
                 fontsize=14)
    ax.tick_params(labelsize=13)
    ax.grid(True, linestyle=":", linewidth=0.8, alpha=0.6)

    class _VegColormapHandler(HandlerBase):
        """Renders a truncated-Greens strip with ticks/labels below it."""
        def create_artists(self, legend, orig_handle,
                           xdescent, ydescent, width, height, fontsize, trans):
            n = 64
            bar_y = ydescent + 0.56 * height
            bar_h = 0.25 * height
            artists = []
            for i in range(n):
                x0 = xdescent + width * i / n
                color = _veg_cmap(_veg_norm(i / (n - 1)))
                artists.append(
                    FancyBboxPatch(
                        (x0, bar_y), width / n + 0.4, bar_h,
                        boxstyle="square,pad=0",
                        linewidth=0, facecolor=color,
                        transform=trans,
                    )
                )
            for frac, lbl in ((0.0, "0"), (0.5, "0.5"), (1.0, "1")):
                xt = xdescent + frac * width
                artists.append(
                    Line2D(
                        [xt, xt],
                        [bar_y - 0.13 * height, bar_y],
                        color="0.2",
                        linewidth=0.8,
                        transform=trans,
                    )
                )
                artists.append(
                    Text(
                        xt,
                        ydescent - 0.03 * height,
                        lbl,
                        ha="center",
                        va="top",
                        fontsize=fontsize * 0.64,
                        color="0.15",
                        transform=trans,
                    )
                )
            return artists

    _handles, _labels = ax.get_legend_handles_labels()
    _handler_map = {}
    if _any_veg_plotted:
        _veg_proxy = Patch(facecolor='none', edgecolor='none')
        _handles.append(_veg_proxy)
        _labels.append('Vegetation cover')
        _handler_map[_veg_proxy] = _VegColormapHandler()

    ax.legend(_handles, _labels, handler_map=_handler_map,
              loc='upper center', bbox_to_anchor=(0.5, -0.18),
              ncol=3, fontsize=11, frameon=True)

    _x_left = 0.0 if x_lim_left is None else float(x_lim_left)
    _y_bottom = 0.0 if y_lim_bottom is None else float(y_lim_bottom)

    ax.set_xlim(left=_x_left, right=x_fill_right)
    if y_lim_top is not None:
        ax.set_ylim(bottom=_y_bottom, top=y_lim_top)
    else:
        ax.set_ylim(bottom=_y_bottom)

    fig.tight_layout()

    plot_filename = f"{output_base}_{timestep}.png"
    fig.savefig(plot_filename, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"Profile snapshot saved to {plot_filename}")
    return plot_filename


def _safe_save_profile_snapshot(dune_x_init, dune_z_init, dune_x_curr, dune_z_curr,
                                substrate_x, substrate_z, output_base, timestep,
                                progress_pct, bc_info=None, vegetation_info=None,
                                y_lim_top=None, y_lim_bottom=None,
                                x_lim_right=None, x_lim_left=None):
    """
    Best-effort wrapper around snapshot plotting that keeps plotting isolated
    from live model state and never interrupts the simulation.
    """
    try:
        return _save_profile_snapshot(
            dict(dune_x_init),
            dict(dune_z_init),
            dict(dune_x_curr),
            dict(dune_z_curr),
            np.asarray(substrate_x, dtype=float).copy(),
            np.asarray(substrate_z, dtype=float).copy(),
            output_base,
            timestep,
            progress_pct,
            bc_info=dict(bc_info) if bc_info is not None else None,
            vegetation_info=dict(vegetation_info) if vegetation_info is not None else None,
            y_lim_top=y_lim_top,
            y_lim_bottom=y_lim_bottom,
            x_lim_right=x_lim_right,
            x_lim_left=x_lim_left,
        )
    except Exception as exc:
        print(
            f"Warning: snapshot plotting failed at timestep {timestep} "
            f"({progress_pct:.1f}%); continuing without this frame. Error: {exc}"
        )
        return None


def _create_evolution_movie(frame_paths, output_path, fps=2):
    """
    Create an animated GIF from a list of snapshot PNG files.

    Parameters
    ----------
    frame_paths : list of str
        Ordered list of PNG file paths to stitch into the animation.
    output_path : str
        Destination path for the GIF file.
    fps : int, optional
        Frames per second (default 2).
    """
    try:
        from PIL import Image
    except ImportError:
        print("Warning: Pillow is not installed – cannot create evolution movie. "
              "Install it with: pip install Pillow")
        return

    images = []
    for path in frame_paths:
        if os.path.isfile(path):
            with Image.open(path) as img:
                images.append(img.copy())

    if not images:
        print("Warning: no snapshot frames found; evolution movie not created.")
        return

    duration_ms = int(1000 / fps)
    images[0].save(
        output_path,
        save_all=True,
        append_images=images[1:],
        duration=duration_ms,
        loop=0
    )
    print(f"Evolution movie saved to {output_path}")


def _create_evolution_video(frame_paths, output_path, fps=2):
    """
    Create an MP4 video from a list of snapshot PNG files.

    Parameters
    ----------
    frame_paths : list of str
        Ordered list of PNG file paths to stitch into the video.
    output_path : str
        Destination path for the MP4 file.
    fps : int, optional
        Frames per second (default 2).
    """
    try:
        import imageio.v2 as imageio
    except ImportError:
        print("Warning: imageio is not installed – cannot create evolution MP4. "
              "Install it with: pip install imageio imageio-ffmpeg")
        return

    valid_paths = [path for path in frame_paths if os.path.isfile(path)]
    if not valid_paths:
        print("Warning: no snapshot frames found; evolution MP4 not created.")
        return

    try:
        with imageio.get_writer(output_path, fps=fps, codec="libx264") as writer:
            for path in valid_paths:
                writer.append_data(imageio.imread(path))
    except Exception as exc:
        print(f"Warning: failed to create evolution MP4 at {output_path}. Error: {exc}")
        return

    print(f"Evolution movie saved to {output_path}")


def run_model(config_path, check_mass_conservation=False,
              plot_evolution=False, plot_evolution_interval=5.0,
              save_evolution_movie=False, skip_first_step_conservation=False,
              manual_plot_limits=False, manual_plot_y_lim_top=None,
              manual_plot_y_lim_bottom=None, manual_plot_x_lim_right=None,
              manual_plot_x_lim_left=None, save_evolution_video=False):
    """
    Run the CS model using a YAML configuration file.

    Parameters
    ----------
    config_path : str
        Path to the YAML configuration file
    check_mass_conservation : bool, optional
        If True, perform mass conservation check after simulation and print report.
        Default is False.
    plot_evolution : bool, optional
        If True, save cross-shore profile evolution snapshots at regular intervals
        during the simulation.  Each snapshot shows the current profile against the
        initial profile (with deposition/erosion shading) and is saved as
        ``<output_base>_<timestep>.png`` in the same directory as the output file.
        Default is False.
    plot_evolution_interval : float, optional
        Interval between snapshots expressed as a percentage of the total simulation
        length (e.g. 5.0 means a snapshot every 5 %).  Only used when
        ``plot_evolution=True``.  Default is 5.0.
    save_evolution_movie : bool, optional
        If True, stitch all saved snapshots into an animated GIF after the
        simulation ends.  The file is saved as ``<output_base>_evolution.gif``
        alongside the model output.  Requires ``plot_evolution=True`` and
        Pillow (``pip install Pillow``).  Default is False.
    save_evolution_video : bool, optional
        If True, stitch all saved snapshots into an MP4 video after the
        simulation ends. The file is saved as ``<output_base>_evolution.mp4``
        alongside the model output. Requires ``plot_evolution=True`` and
        imageio/imageio-ffmpeg (``pip install imageio imageio-ffmpeg``).
        Default is False.
    skip_first_step_conservation : bool, optional
        If True, the first model time step is excluded from the mass conservation
        analysis.  Passed directly to ``check_mass_conservation`` when
        ``check_mass_conservation=True``.  Default is False.
    manual_plot_limits : bool, optional
        If True, manually provided plot limits are used for evolution snapshots.
        If False, automatic limits are used.  Default is False.
    manual_plot_y_lim_top : float or None, optional
        Manual upper y-axis limit (meters) for evolution snapshots when
        ``manual_plot_limits=True``. If None, automatic y-limit is kept.
        Default is None.
    manual_plot_y_lim_bottom : float or None, optional
        Manual lower y-axis limit (meters) for evolution snapshots when
        ``manual_plot_limits=True``. If None, defaults to 0.
        Default is None.
    manual_plot_x_lim_right : float or None, optional
        Manual right x-axis limit (meters) for evolution snapshots when
        ``manual_plot_limits=True``. If None, automatic x-limit is kept.
        Default is None.
    manual_plot_x_lim_left : float or None, optional
        Manual left x-axis limit (meters) for evolution snapshots when
        ``manual_plot_limits=True``. If None, defaults to 0.
        Default is None.

    Returns
    -------
    results : dict
        Model output
    """

    # Load configuration
    config = load_config(config_path)


    # Inject each config key as a variable in the global namespace
    # for section, params in config.items():
    #     if isinstance(params, dict):
    #         for k, v in params.items():
    #             globals()[k] = v
    #     else:
    #         globals()[section] = params

    params = {}
    for section, values in config.items():
        if isinstance(values, dict):
            params.update(values)
        else:
            params[section] = values

    p = SimpleNamespace(**params)

    # load boundary conditions
    time_length_dt = int(p.time_length_hours * 3600 / p.dt)

    wave_time_series = load_timeseries(
        p.wave_height, time_length_dt, column="wave_height"
    )

    wave_period_time_series = load_timeseries(
        p.wave_period, time_length_dt, column="wave_period"
    )

    wave_dir_time_series = load_timeseries(
        p.wave_direction, time_length_dt, column="wave_direction"
    )

    swl_time_series = load_timeseries(
        p.swl, time_length_dt, column="swl"
    )

    wind_speed_time_series = load_timeseries(
        p.wind_speed, time_length_dt, column="wind_speed"
    )

    wind_dir_time_series = load_timeseries(
        p.wind_direction, time_length_dt, column="wind_direction"
    )

    month_timeseries = load_timeseries(
        p.month, time_length_dt, column="month", default=7
    )


    # Allocate memory
    qbar = 0#np.zeros(time_length_hours)
    qnour = 0#np.zeros(time_length_hours) #onshore transport from shoreface nourishments
    qslr = 0#np.zeros(time_length_hours)
    dvol = 0#np.zeros(time_length_hours)

    IV = 1 #counter for nourishments
    vnourd = np.zeros(time_length_dt) #nourished volume dune
    vnourb = np.zeros(time_length_dt) #nourished volume beach
    vnourbar = np.zeros(time_length_dt) #nourished volume bar

    volat = np.asarray(getattr(p, 'volat', []), dtype=int)   # time step when nourishment event k+1 happens
    volab = np.asarray(getattr(p, 'volab', []), dtype=float)  # beach volume added at event k+1
    volad = np.asarray(getattr(p, 'volad', []), dtype=float)  # dune volume added at event k+1
    volrb = np.asarray(getattr(p, 'volrb', []), dtype=float)  # bar volume added at event k+1
    
    vanr_discrete = len(volat)
    beach_rate = getattr(p, "beach_nourishment_rate", 0.0)
    
    if beach_rate > 0.0:
        # Continuous nourishment: one event per timestep
        volat_cont = np.arange(time_length_dt, dtype=int)
        volab_cont = np.full(time_length_dt, beach_rate, dtype=float)
        volad_cont = np.zeros(time_length_dt, dtype=float)
        volrb_cont = np.zeros(time_length_dt, dtype=float)
    
        # Add discrete nourishment on top
        for k in range(vanr_discrete):
            t = volat[k]
            if 0 <= t < time_length_dt:
                volab_cont[t] += volab[k]
                volad_cont[t] += volad[k]
                volrb_cont[t] += volrb[k]

        # Replace nourishment schedule
        volat = volat_cont
        volab = volab_cont
        volad = volad_cont
        volrb = volrb_cont
        p.vanr = time_length_dt
        IV = 1
    
    barnour = np.ones(time_length_dt)*p.barnour_total
    eroded_volume_time_series = np.zeros(time_length_dt)
    overwash_volume_time_series = np.zeros(time_length_dt)
    dune_height_time_series = np.zeros(time_length_dt)
    dune_foot_time_series = np.zeros(time_length_dt)
    total_water_level_time_series = np.zeros(time_length_dt)
    runup_time_series = np.zeros(time_length_dt)
    X_time_series = np.zeros(time_length_dt)
    vbeach_volume_time_series = np.zeros(time_length_dt)

    ##### Define initial profile #####

    # define non-erodible substrate (dike + flat beach base)
    substrate_x = [0, p.yl_hs, p.ylprim_hs, p.ysprim_hs, p.ys_hs, p.yg]
    substrate_z = [p.base_elevation_back, p.base_elevation_back, p.s_hs + p.base_elevation_front, p.s_hs + p.base_elevation_front, p.base_elevation_front, p.base_elevation_front] # use base 0 = mean sea level throughout the profile for the calculations

    # Define dune profile as dictionaries for clarity
    dune_x = {
        'back_ref': 0.0,
        'lee_foot': p.yl,
        'lee_crest': p.ylprim,
        'crest_temp': p.ysprim,
        'seaward_crest': p.ysprim,
        'ramp_dike_temp': p.ysprim, # temporary point to help with calculations
        'ramp_crest': p.ysprim,  # initially same as seaward_crest
        'ramp_toe_temp': p.ys,
        'ramp_toe': p.ys,
        'shoreline': p.yg
    }

    lee_base = getattr(p, "dfoot_back", p.dfoot)
    dune_z = {
        'back_ref': lee_base,
        'lee_foot': lee_base,
        'lee_crest': p.s + p.dfoot,
        'crest_temp': p.s + p.dfoot,
        'seaward_crest': p.s + p.dfoot,
        'ramp_dike_temp': p.sramp + p.dfoot, # temporary point to help with calculations
        'ramp_crest': p.sramp + p.dfoot,
        'ramp_toe_temp': p.dfoot,
        'ramp_toe': p.dfoot,
        'shoreline': 0.0
    }
    
    # Initialize results dictionary to store time series
    results = {
        'dune_x': [dune_x.copy()],
        'dune_z': [dune_z.copy()],
        'substrate_x': substrate_x,
        'substrate_z': substrate_z,
        'vramp': [],
        'vbeach': [],
        'vdune': [],
        'wave_height': [],
        'swl': [],
        'runup': [],
        'total_water_level': [],
        'eroded_volume': [],
        'overwash_volume': [],
        'dune_height': [dune_z['seaward_crest'] - dune_z['ramp_toe']],
        'dune_foot': [dune_x['ramp_toe']],
        'rhov_df': [],
        'rhov_dt': [],
        'rhov_db': [],
        'rhov_df2': [],
        'rhov_dt2': [],
        'rhov_db2': [],
        'volbar': [],
        'qs_crest': [],
        'qs_lee': [],
        'qs_ramp': [],
        'qwnet': [],
        'qbar': [],
        'wave_period': wave_period_time_series,
        'wave_direction': wave_dir_time_series,
        'wind_speed': wind_speed_time_series,
        'wind_direction': wind_dir_time_series,
        'time': np.arange(0, time_length_dt*p.dt, p.dt),
        # Conservation time-series
        'conservation': {
            'dune_error': [],
            'beach_error': [],
            'bar_error': [],
            'total_error': [],
            'dune_error_pct': [],
            'beach_error_pct': [],
            'bar_error_pct': [],
            'total_error_pct': [],
            'dune_volume': [],
            'beach_volume': [],
            'bar_volume': [],
            'total_volume': []
        },
        # Nourishment volumes applied each timestep (external sediment sources)
        'nourishment': {
            'dune': [],
            'beach': [],
            'bar': []
        }
    }

    # Initialize cumulative time for vegetation update
    cumdt = 0
    vegupdate = 0
    s_veg = p.s # allocate value for the vegetation growth calc
    ys_veg = p.ys # allocate value for the vegetation growth calc
    yl_veg = p.yl # allocate value for the vegetation growth calc
    sramp_veg = p.sramp # allocate value for the vegetation growth calc

    # Initialize varying parameters
    # Calculate initial vbeach from shoreline position (yg) instead of using config value
    # to ensure consistency with the beach volume formula
    vbeach = ((p.yg - p.ys) - p.bl) / p.al
    vramp = p.vramp
    rhov_df = p.rhov_df
    rhov_dt = p.rhov_dt
    rhov_db = p.rhov_db
    rhov_df2 = p.rhov_df2
    rhov_dt2 = p.rhov_dt2
    rhov_db2 = p.rhov_db2
    volbar = p.volbar
    deltavall = 0

    # Calculate initial volumes for conservation tracking
    initial_volumes = calculate_compartment_volumes(
        dune_x, dune_z, substrate_x, substrate_z, vbeach, volbar, p.dclos, p.dfoot
    )
    
    # Store initial volumes in conservation tracking
    for key in ['dune_volume', 'beach_volume', 'bar_volume', 'total_volume']:
        compartment = key.replace('_volume', '')
        results['conservation'][key].append(initial_volumes[compartment])
    
    # Initialize with zero errors for initial timestep
    for key in ['dune_error', 'beach_error', 'bar_error', 'total_error',
                'dune_error_pct', 'beach_error_pct', 'bar_error_pct', 'total_error_pct']:
        results['conservation'][key].append(0.0)

    # Initialize nourishment series with zero for the initial (pre-loop) timestep
    for key in ['dune', 'beach', 'bar']:
        results['nourishment'][key].append(0.0)

    # Store reference to previous volumes for conservation tracking
    prev_volumes = initial_volumes.copy()

    # Set up evolution plotting (resolve output base early so snapshots can be saved
    # incrementally during the loop)
    output_filename = config.get('output', {}).get('filename', 'cs_model_results.pkl')
    if plot_evolution:
        _plot_output_base = os.path.splitext(output_filename)[0]
        _plot_next_pct = float(plot_evolution_interval)
        _total_steps = len(wave_time_series)
        _plot_saved_frames = []
        # Erosion-event snapshot tracking
        _erosion_active = False
        _erosion_event_steps = 0
        _erosion_mid_saved = False

        # Pre-compute a consistent y-axis upper limit from the full boundary
        # time series so that snapshots do not rescale when large water-level
        # events occur.  Uses the initial profile slope as a representative
        # estimate of runup.
        _slope_init = dune_z['ramp_toe'] / max(
            dune_x['shoreline'] - dune_x['ramp_toe'], 1e-10
        )
        _h0_all  = np.array(wave_time_series) / 1.414  # Hs -> H_rms (H_rms = Hs/sqrt(2))
        _tp_all  = np.array(wave_period_time_series)
        _swl_all = np.array(swl_time_series)
        # Deep-water wavelength: L0 = g*T^2/(2*pi) ≈ 1.5613*(0.95*Tp)^2
        _LO_all  = 1.5613 * (0.95 * _tp_all) ** 2
        # Runup estimate (Hallin et al. 2019): r = (0.158/0.03)*slope*sqrt(H*L0)
        _r_all   = (0.158 / 0.03) * _slope_init * np.sqrt(_h0_all * _LO_all)
        _twl_all = _r_all + _swl_all
        # Add 15% headroom above the peak TWL for visual clarity
        _plot_y_max = float(np.nanmax(_twl_all)) * 1.15
        # Also ensure the initial dune crest (+ 0.5 m buffer) fits within the axis
        _plot_y_max = max(_plot_y_max, (dune_z['seaward_crest'] + 6.0) * 1.15)
        # Pre-compute a fixed x-axis right limit from the initial profile so
       # that all GIF frames share the same horizontal extent.  Use 20 %
        # margin beyond the initial shoreline (same shape as the dynamic formula but evaluated only once).
        _plot_x_max = dune_x['shoreline'] + 0.2 * (dune_x['shoreline'] - dune_x['ramp_toe'])
        _plot_y_min = 0.0
        _plot_x_min = 0.0

        if manual_plot_limits:
            if manual_plot_y_lim_top is not None:
                if float(manual_plot_y_lim_top) <= 0.0:
                    raise ValueError("manual_plot_y_lim_top must be > 0 when provided.")
                _plot_y_max = float(manual_plot_y_lim_top)
            if manual_plot_y_lim_bottom is not None:
                _plot_y_min = float(manual_plot_y_lim_bottom)
            if manual_plot_x_lim_right is not None:
                if float(manual_plot_x_lim_right) <= 0.0:
                    raise ValueError("manual_plot_x_lim_right must be > 0 when provided.")
                _plot_x_max = float(manual_plot_x_lim_right)
            if manual_plot_x_lim_left is not None:
                _plot_x_min = float(manual_plot_x_lim_left)

            if _plot_y_max <= _plot_y_min:
                raise ValueError(
                    "manual y-axis limits are invalid: manual_plot_y_lim_top "
                    "must be greater than manual_plot_y_lim_bottom."
                )
            if _plot_x_max <= _plot_x_min:
                raise ValueError(
                    "manual x-axis limits are invalid: manual_plot_x_lim_right "
                    "must be greater than manual_plot_x_lim_left."
                )
        elif (
            (manual_plot_y_lim_top is not None)
            or (manual_plot_y_lim_bottom is not None)
            or (manual_plot_x_lim_right is not None)
            or (manual_plot_x_lim_left is not None)
        ):
            print(
                "Warning: manual plot limits were provided "
                "but manual_plot_limits=False; automatic plot limits are used."
            )
        
    ############################################################################
    # START SIMULATION
    ############################################################################

    for i in range(len(wave_time_series)):

        month_i = month_timeseries[i]
        ### update counter for vegetation
        cumdt = cumdt + p.dt
        if cumdt >= p.vegdt:
            # Update vegetation parameters here
            cumdt = 0
            vegupdate = 1

        ### hydrodynamics ###
        h0 = wave_time_series[i] / 1.414 #rms wave height
        tp = wave_period_time_series[i]  # Wave period (s)
        swl = swl_time_series[i]  # Sea water level (m)

        # modify wave height with respect to wave direction
        if 90 <= p.shorn <= 270:
            UL = p.shorn + 90
            LL = p.shorn - 90
            if LL <= wave_dir_time_series[i] <= UL:
                teta = wave_dir_time_series[i] - p.shorn
            else:
                teta = 0.0
                h0 = 0.0
                t = 0.0
        elif p.shorn > 270:
            UL = p.shorn + 90 - 360
            LL = p.shorn - 90
            if wave_dir_time_series[i] >= LL:
                teta = wave_dir_time_series[i] - p.shorn
            elif wave_dir_time_series[i] <= UL:
                teta = wave_dir_time_series[i] + 360 - p.shorn
            else:
                teta = 0.0
                h0 = 0.0
                t = 0.0
        elif p.shorn < 90:
            UL = p.shorn + 90
            LL = p.shorn - 90 + 360
            if wave_dir_time_series[i] <= UL:
                teta = wave_dir_time_series[i] - p.shorn
            elif wave_dir_time_series[i] >= LL:
                teta = 360 - wave_dir_time_series[i] + p.shorn
            else:
                teta = 0.0
                h0 = 0.0
                tp = 0.0

        # Do not allow the waves to become too small
        if tp <= 0.0:
            tp = 1.0
        if h0 <= 0.0:
            h0 = 0.1

        # calculate approximate runup
        
        LO = 1.5613 * (0.95 * tp) ** 2
        slope_foreshore = dune_z['ramp_toe']/(dune_x['shoreline']- dune_x['ramp_toe'])
        r = (0.158/0.03)*slope_foreshore*np.sqrt(h0 * LO) #we divide by an actual slope and recalculate based on the slope we have so to vary potential dissipation with X
        X = 2 * vbeach / p.dfoot * (1 - swl / p.dfoot)

        if X < 0:
            X = 0 # surge above the dunefoot, hence no bed for runup friction, only wave dissipation.
            
        rprim = (r) * np.exp(-2 * p.cf * X) + (p.dfoot) * (1 - np.exp(-2 * p.cf * X))
            
        if rprim > p.dfoot:
            r = rprim

        runup_time_series[i] = r
        X_time_series[i] = X
        
        aover = 3.0
        alfa = 0.0
        # calculate potential erosion

        # total frontal area per alongshore cm (cm^2/cm)
        Av_cm2_per_cm_spec1 = max(0.0, p.frontal_cover_scale_spec1*rhov_df + p.frontal_cover_scale_spec1*rhov_dt + p.frontal_cover_scale_spec1*rhov_db)
        Av_cm2_per_cm_spec2 = max(0.0, p.frontal_cover_scale_spec2*rhov_df2 + p.frontal_cover_scale_spec2*rhov_dt2 + p.frontal_cover_scale_spec2*rhov_db2)

        Av_total = Av_cm2_per_cm_spec1 + Av_cm2_per_cm_spec2
        # Figlus-like slope: 2.1% reduction per (cm^2/cm) -> 0.021 fraction per (cm^2/cm)
        # multiplier to apply to qd (1 = no effect, 0.72 = 28% reduction)
        vegred = (1.0 - min(0.021 * Av_total, 0.28))

        if r + swl > p.dfoot:
            if dune_z['seaward_crest'] >= dune_z['ramp_dike_temp']:
                if r + swl <= dune_z['seaward_crest']:
                    qd = vegred * 4 * p.cimpact * (r + swl - p.dfoot) ** 2 / tp
                else:
                    print('OVERWASH')
                    qd = vegred * 4 * p.cimpact * (((r + swl - p.dfoot) * (dune_z['seaward_crest']-p.dfoot)) / tp)
                    alfa = (((r + swl - p.dfoot) / (dune_z['seaward_crest']-p.dfoot)) - 1) / aover
            else: 
                if r + swl <= dune_z['ramp_dike_temp']:
                    qd = vegred * 4 * p.cimpact * (r + swl - p.dfoot) ** 2 / tp
                else:
                    print('OVERWASH')
                    qd = vegred * 4 * p.cimpact * (((r + swl - p.dfoot) * (dune_z['ramp_dike_temp']-p.dfoot)) / tp)
                    alfa = (((r + swl - p.dfoot) / (dune_z['ramp_dike_temp']-p.dfoot)) - 1) / aover
        else:
            qd = 0.0
        
        qs = qd / (1 + alfa)
        ql = qd * alfa / (1 + alfa)

        total_water_level = r + swl 

        total_water_level_time_series[i] = total_water_level

        #### aeolian transport at the beach ####
            ## calculate wind direction to the beach normal
        if 80 <= p.shorn <= 280:
            UL = p.shorn + 80
            LL = p.shorn - 80
            if LL <= wind_dir_time_series[i] <= UL:
                tetaw = wind_dir_time_series[i] - p.shorn
                wspeed = wind_speed_time_series[i]
            else:
                tetaw = 0.0
                wspeed = 0.0
        elif p.shorn > 280:
            UL = p.shorn + 80 - 360
            LL = p.shorn - 80
            if wind_dir_time_series[i] >= LL:
                tetaw = wind_dir_time_series[i] - p.shorn
                wspeed = wind_speed_time_series[i]
            elif wind_dir_time_series[i] <= UL:
                tetaw = wind_dir_time_series[i] + 360 - p.shorn
                wspeed = wind_speed_time_series[i]
            else:
                tetaw = 0.0
                wspeed = 0.0
        elif p.shorn < 80:
            UL = p.shorn + 80
            LL = p.shorn - 80 + 360
            if wind_dir_time_series[i] <= UL:
                tetaw = wind_dir_time_series[i] - p.shorn
                wspeed = wind_speed_time_series[i]
            elif wind_dir_time_series[i] >= LL:
                tetaw = 360 - wind_dir_time_series[i] + p.shorn
                wspeed = wind_speed_time_series[i]
            else:
                tetaw = 0.0
                wspeed = 0.0
        
        wdir = tetaw # final adjusted wind direction

        ## Calculate wind shear and critical shear at the beach
        wcrit = p.awind * np.sqrt((p.rhos - p.rhoa) / p.rhoa * p.g * p.d50)
        #z0 = 0.081 * np.log10(d50 * 10 ** 3 / 0.18) # the one used in Hallin's model
        z0 = p.d50/30 # the one used in Bagnold

        #wshear = wspeed * max(kappa / np.log(z / z0), 0.053)
        wshear = wspeed * p.kappa / np.log(p.z / z0)  # shear velocity

        if wshear > wcrit:
            mwe = p.kw * np.sqrt(p.d50 / p.d50ref) * p.rhoa * wshear ** 2 / p.g * (wshear - wcrit)
        else:
            mwe = 0
        qwe = mwe / (p.rhos * (1 - p.p))

        #yr = ys + (1 - (r + swl) / dfoot) * (yg - ys)
        yr = dune_x['ramp_toe'] + (1 - (r + swl) / p.dfoot) * (dune_x['shoreline'] - dune_x['ramp_toe'])
        #bdry = yr - ys
        bdry = yr - dune_x['ramp_toe']

        if r + swl >= p.dfoot or yr < dune_x['ramp_toe']:
            bdry = 0
        f = bdry / abs(np.cos(np.radians(wdir)))
        qw = qwe * (1 - np.exp(-p.cwrate * f))
        qwnet = qw * np.cos(np.radians(wdir))

        # Introduce dynamic reduction based on Duran and Moore (2013) Ad
        #qwnet = qwnet * max(1.0 - 4.2 * (s-p.dfoot) / (yg - ys), 0.0)
        if dune_z['seaward_crest'] >= dune_z['ramp_dike_temp']:
            qwnet = qwnet * max(1.0 - 4.2 * (dune_z['seaward_crest']-p.dfoot) / (dune_x['shoreline'] - dune_x['ramp_toe']), 0.0)
        else:
            qwnet = qwnet * max(1.0 - 4.2 * (dune_z['ramp_dike_temp']-p.dfoot) / (dune_x['shoreline'] - dune_x['ramp_toe']), 0.0)

        ######## AEOLIAN DEPOSITION AT DUNE ########

        #if i>=4464:
        #    print(i)
            
        qs_ramp, qs_crest, qs_lee = deposition_dune(dune_x, dune_z, substrate_x, substrate_z, qwnet, p.dt, tetaw, wspeed, p.rhos, p.rhoa, p.d50, p.kappa, p.z, p.vegcount, p.awind, rhov_df=rhov_df, rhov_dt=rhov_dt, rhov_db=rhov_db, 
            rhov_df2=rhov_df2 if p.vegcount == 1 else 0, 
            rhov_dt2=rhov_dt2 if p.vegcount == 1 else 0, 
            rhov_db2=rhov_db2 if p.vegcount == 1 else 0)

        ######## NOURISHMENT AND SUBTIDAL EVOLUTION ########
        
        wfall = FALVEL(p.temp, p.d50) # calculate fall velocity
        
        # subtidal evolution
        qbar, qnour, qslr, dvol = subaq(p.dt, volbar, h0, tp, wfall, p.slr, barnour[i], p.cb, p.bactive)

        # specify nourished volumes
        if p.vanr > 0:
            vnourd[i], vnourb[i], vnourbar[i], IV = nour(i, IV, volat, volad, volab, volrb, p.vanr)
            
        # Calculate bar evolution
        volbar, barnour[i] = bar(volbar, qbar, p.dt, vnourbar[i], barnour[i], qnour)
        
        # Calculate beach volume BEFORE morpho_update (using old ramp_toe)
        vbeach_old = ((dune_x['shoreline'] - dune_x['ramp_toe'])-p.bl)/p.al
        db_old = vbeach_old / (dune_x['shoreline'] - dune_x['ramp_toe'])
        vbeachtot_old = dune_x['ramp_toe']  * (p.dclos + p.dfoot) + (dune_x['shoreline'] - dune_x['ramp_toe']) * (p.dclos + db_old)
        
        if swl>0:
            print(i)
            
        ####### MORPHOLOGICAL UPDATE OF THE DUNE PROFILE #########
        dune_x, dune_z, deltavramp, erod_overwash, deltavramp_ero, deltavall = morpho_update(dune_x, dune_z, substrate_x, substrate_z, qs_ramp, qs_crest, qs_lee, qs, ql, r, swl, vnourd[i], p.dt, p.betas, p.betas_hs, p.betal, p.betal_hs, vramp, deltavall)

        ######## VEGETATION GROWTH UPDATE #########
            # update vegetation parameters
        if p.vegcount == 1 and vegupdate == 1:
            rhov_df, rhov_dt, rhov_db, rhov_df2, rhov_dt2, rhov_db2 = vegNield_seasonal(p.vegdt,s_veg,dune_z['seaward_crest'], dune_x['lee_foot'], yl_veg, dune_x['ramp_toe'], ys_veg, dune_z['ramp_crest'], sramp_veg, p.betas, p.betal, rhov_df, rhov_dt, rhov_db, rhov_df2, rhov_dt2, rhov_db2, month_i)
            vegupdate = 0
            s_veg = dune_z['seaward_crest']
            ys_veg = dune_x['ramp_toe']
            yl_veg = dune_x['lee_foot']
            sramp_veg = dune_z['ramp_crest']
        elif p.vegcount == 0 and vegupdate == 1:
            rhov_df, rhov_dt, rhov_db = veg(ys, ysnew, yl, ylnew, sramp, srampnew, betas, betal, s, snew, vegdt,rhov_df,rhov_dt,rhov_db)
            vegupdate = 0

        
        # update values for next iteration
        vramp = partial_volume(dune_x, dune_z, substrate_x, substrate_z, dune_x['ramp_crest'], dune_x['ramp_toe'])
        vdune = partial_volume(dune_x, dune_z, substrate_x, substrate_z, dune_x['lee_foot'], dune_x['ramp_toe'])
        dune_height_time_series[i] = dune_z['seaward_crest'] - dune_z['ramp_toe']
        dune_foot_time_series[i] = dune_x['ramp_toe']
        total_water_level_time_series[i] = total_water_level
        eroded_volume_time_series[i] = deltavramp_ero #erod_ramp
        overwash_volume_time_series[i] = erod_overwash
                
        # Calculate beach volume after morpho_update, but before the shoreline movement (so to not account the volume created by the moving dune)
        vbeach_phantom = ((dune_x['shoreline'] - dune_x['ramp_toe'])-p.bl)/p.al
        db_phantom = vbeach_phantom / (dune_x['shoreline'] - dune_x['ramp_toe'])
        vbeachtot_phantom = dune_x['ramp_toe']  * (p.dclos + p.dfoot) + (dune_x['shoreline'] - dune_x['ramp_toe']) * (p.dclos + db_old)
        
        # Calculate new beach volume including eroded sediment from dune
        vbeachtotnew = vbeachtot_phantom + (-qwnet - qslr + qnour - qbar) * p.dt + vnourb[i] + deltavramp_ero
        vbeachtotmin = dune_x['ramp_toe'] * (p.dclos + p.dfoot) + 0.5 * p.dfoot**2 / p.betaf + p.dfoot / p.betaf * p.dclos

        # Calculate new shoreline coordinate
        dune_x['shoreline'] = dune_x['ramp_toe'] + (p.al * (vbeachtotnew - dune_x['ramp_toe'] * (p.dfoot + p.dclos)) + p.bl) / (1 + p.al * p.dclos)
        
        # Calculate final beach volume for output
        vbeach = ((dune_x['shoreline'] - dune_x['ramp_toe'])-p.bl)/p.al
        vbeach_volume_time_series[i] = vbeach # actual vbeach
        
        #delta_vbeach_volume_time_series[i] = vbeach - vbeach_phantom #change based on the phantom
        #delta_vbeachtot_volume_time_series[i] = vbeachtotnew - vbeachtot_phantom #beachtot change based on phantom

        ####### CALCULATE CONSERVATION METRICS #######
        # Calculate current volumes
        current_volumes = calculate_compartment_volumes(
            dune_x, dune_z, substrate_x, substrate_z, vbeach, volbar, p.dclos, p.dfoot
        )
        
        # Prepare flux dictionary for conservation calculation
        fluxes = {
            'qs_ramp': qs_ramp,
            'qs_crest': qs_crest,
            'qs_lee': qs_lee,
            'qwnet': qwnet,
            'qbar': qbar,
            'eroded_volume': deltavramp_ero,
            'nour_dune': vnourd[i],
            'nour_beach': vnourb[i],
            'nour_bar': vnourbar[i]
        }

        # Store nourishment volumes applied this timestep
        results['nourishment']['dune'].append(vnourd[i])
        results['nourishment']['beach'].append(vnourb[i])
        results['nourishment']['bar'].append(vnourbar[i])
        
        # Calculate conservation metrics for this timestep
        conservation_metrics = calculate_timestep_conservation(
            current_volumes, prev_volumes, fluxes, p.dt, initial_volumes
        )
        
        # Store conservation metrics
        for key, value in conservation_metrics.items():
            results['conservation'][key].append(value)
        
        # Update previous volumes for next iteration
        prev_volumes = current_volumes.copy()
        
        ####### STORE OUTPUTS #######
        results['dune_x'].append(dune_x.copy())
        results['dune_z'].append(dune_z.copy())
        results['vramp'].append(vramp)
        results['vbeach'].append(vbeach)
        results['vdune'].append(vdune)
        results['wave_height'].append(wave_time_series[i])
        results['swl'].append(swl)
        results['runup'].append(r)
        results['total_water_level'].append(total_water_level)
        results['eroded_volume'].append(deltavramp_ero)
        results['overwash_volume'].append(erod_overwash)
        results['dune_height'].append(dune_z['seaward_crest'] - dune_z['ramp_toe'])
        results['dune_foot'].append(dune_x['ramp_toe'])
        results['rhov_db'].append(rhov_db)
        results['rhov_dt'].append(rhov_dt)
        results['rhov_df'].append(rhov_df)
        results['rhov_db2'].append(rhov_db2)
        results['rhov_dt2'].append(rhov_dt2)
        results['rhov_df2'].append(rhov_df2)
        results['volbar'].append(volbar)
        results['qs_crest'].append(qs_crest)
        results['qs_lee'].append(qs_lee)
        results['qs_ramp'].append(qs_ramp)
        results['qwnet'].append(qwnet)
        results['qbar'].append(qbar)

        ## print progress
        if i % 1000 == 0:
            print(f"Progress: {100*i/len(wave_time_series):.1f}% ({i} of {len(wave_time_series)} time steps)")

        ## save profile evolution snapshot (if requested)
        if plot_evolution:
            current_pct = 100.0 * (i + 1) / _total_steps

            # Build shared bc_info for all snapshot calls this timestep
            _bc_info = {
                'wave_height': wave_time_series[i],
                'wave_period': tp,
                'swl': swl,
                'runup': r,
                'total_water_level': total_water_level,
                'wind_speed': wind_speed_time_series[i],
                'wind_direction': wind_dir_time_series[i],
            }
            _veg_info = {
                'rhov_df': rhov_df,
                'rhov_dt': rhov_dt,
                'rhov_db': rhov_db,
                'rhov_df2': rhov_df2 if p.vegcount == 1 else 0.0,
                'rhov_dt2': rhov_dt2 if p.vegcount == 1 else 0.0,
                'rhov_db2': rhov_db2 if p.vegcount == 1 else 0.0,
            }

            # ---- Erosion-event snapshots ----
            _is_erosive = (deltavramp_ero > 0)
            if _is_erosive and not _erosion_active:
                # First timestep of a new erosion event
                _erosion_active = True
                _erosion_event_steps = 1
                _erosion_mid_saved = False
                _fp = _safe_save_profile_snapshot(
                    results['dune_x'][0], results['dune_z'][0],
                    dune_x, dune_z,
                    substrate_x, substrate_z,
                    _plot_output_base, i + 1, current_pct,
                    bc_info=_bc_info, vegetation_info=_veg_info,
                    y_lim_top=_plot_y_max, y_lim_bottom=_plot_y_min,
                    x_lim_right=_plot_x_max, x_lim_left=_plot_x_min
                )
                if _fp is not None:
                    _plot_saved_frames.append(_fp)
            elif _is_erosive and _erosion_active:
                _erosion_event_steps += 1
                if not _erosion_mid_saved and _erosion_event_steps >= 10:
                    # Mid-event snapshot (after at least 10 erosive steps)
                    _erosion_mid_saved = True
                    _fp = _safe_save_profile_snapshot(
                        results['dune_x'][0], results['dune_z'][0],
                        dune_x, dune_z,
                        substrate_x, substrate_z,
                        _plot_output_base, i + 1, current_pct,
                        bc_info=_bc_info, vegetation_info=_veg_info,
                        y_lim_top=_plot_y_max, y_lim_bottom=_plot_y_min,
                        x_lim_right=_plot_x_max, x_lim_left=_plot_x_min
                    )
                    if _fp is not None:
                        _plot_saved_frames.append(_fp)
            elif not _is_erosive and _erosion_active:
                # Last timestep of the erosion event – save end-of-event snapshot
                _erosion_active = False
                _fp = _safe_save_profile_snapshot(
                    results['dune_x'][0], results['dune_z'][0],
                    dune_x, dune_z,
                    substrate_x, substrate_z,
                    _plot_output_base, i + 1, current_pct,
                    bc_info=_bc_info, vegetation_info=_veg_info,
                    y_lim_top=_plot_y_max, y_lim_bottom=_plot_y_min,
                    x_lim_right=_plot_x_max, x_lim_left=_plot_x_min
                )
                if _fp is not None:
                    _plot_saved_frames.append(_fp)

            # ---- Regular interval snapshots ----
            if current_pct >= _plot_next_pct:
                _frame_path = _safe_save_profile_snapshot(
                    results['dune_x'][0], results['dune_z'][0],
                    dune_x, dune_z,
                    substrate_x, substrate_z,
                    _plot_output_base, i + 1, current_pct,
                    bc_info=_bc_info, vegetation_info=_veg_info,
                    y_lim_top=_plot_y_max, y_lim_bottom=_plot_y_min,
                    x_lim_right=_plot_x_max, x_lim_left=_plot_x_min
                )
                if _frame_path is not None:
                    _plot_saved_frames.append(_frame_path)
                _plot_next_pct += float(plot_evolution_interval)

    ############################################################################
    # SAVE RESULTS
    ############################################################################

    with open(output_filename, 'wb') as f:
        pickle.dump(results, f)

    print(f"Simulation completed. Results saved to {output_filename}.")

    ############################################################################
    # EVOLUTION MOVIE (if requested)
    ############################################################################

    if save_evolution_movie and not plot_evolution:
        print("Warning: save_evolution_movie=True has no effect when plot_evolution=False.")
    if save_evolution_video and not plot_evolution:
        print("Warning: save_evolution_video=True has no effect when plot_evolution=False.")

    if plot_evolution and save_evolution_movie and _plot_saved_frames:
        movie_path = f"{_plot_output_base}_evolution.gif"
        try:
            _create_evolution_movie(_plot_saved_frames, movie_path)
        except Exception as exc:
            print(f"Warning: evolution GIF creation failed; continuing without GIF. Error: {exc}")
    if plot_evolution and save_evolution_video and _plot_saved_frames:
        video_path = f"{_plot_output_base}_evolution.mp4"
        try:
            _create_evolution_video(_plot_saved_frames, video_path)
        except Exception as exc:
            print(f"Warning: evolution MP4 creation failed; continuing without MP4. Error: {exc}")
    
    ############################################################################
    # MASS CONSERVATION CHECK (if requested)
    ############################################################################
    
    if check_mass_conservation:
        print()  # Empty line for better formatting
        from cs_model.mass_conservation import check_mass_conservation as check_mass
        mass_check_result = check_mass(results, config, verbose=True,
                                       skip_first_step=skip_first_step_conservation)
        
        # Save the full mass conservation report to file
        if mass_check_result:
            mass_report_filename = output_filename.replace('.pkl', '_mass_conservation.txt')
            with open(mass_report_filename, 'w') as f:
                f.write("Detailed Mass Conservation Report\n")
                f.write("=" * 80 + "\n")
                f.write(f"Output file: {output_filename}\n")
                f.write("=" * 80 + "\n")
                f.write("\n")
                # Write the full detailed report
                f.write(mass_check_result['report_text'])
            print(f"Mass conservation report saved to {mass_report_filename}")

            # Print final summary to terminal for visibility
            print()
            print("=" * 80)
            print("MASS CONSERVATION SUMMARY")
            print("=" * 80)
            print(f"Status: {'✓ PASSED' if mass_check_result['conserved'] else '⚠ FAILED'}")
            print(f"Total Volume Error:         {mass_check_result['total_error']:+.6f} m³/m ({mass_check_result['total_error_percent']:.4f}% of V₀)")
            print(f"Transport-relative error:   {100*mass_check_result['transport_relative_error']:.4f}%")
            print(f"Annualised drift rate:       {mass_check_result['drift_per_year']:+.3e} m³/m/year")
            print(f"Initial Total Volume:        {mass_check_result['initial_total']:.4f} m³/m")
            print(f"Final Total Volume:          {mass_check_result['final_total']:.4f} m³/m")
            print(f"Thresholds — transport: {mass_check_result['transport_threshold']}% | volume: {mass_check_result['volume_threshold']}%")
            if mass_check_result['conserved']:
                print("\n✓ Mass is conserved within acceptable limits!")
            else:
                print("\n⚠ WARNING: Mass balance error exceeds acceptable threshold!")
            print("=" * 80)
            print()

    return results

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Usage: python main.py path/to/config.yaml")

    run_model(sys.argv[1])
