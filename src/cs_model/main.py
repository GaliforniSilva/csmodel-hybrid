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
from shapely.geometry import LineString
from scipy.interpolate import interp1d
import pandas as pd
import netCDF4 as nc
import yaml
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
                           progress_pct, bc_info=None):
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
        (all optional): ``swl`` (m), ``wave_height`` (Hs, m),
        ``wave_period`` (Tp, s), ``runup`` (m), ``total_water_level`` (m).
        When ``swl`` is present a water-level line is drawn on the seaward
        side of the profile.

    Returns
    -------
    plot_filename : str or None
        Path of the saved PNG file, or ``None`` if the snapshot was skipped.
    """
    c_init  = "#1f77b4"   # blue
    c_final = "#d62728"   # red
    c_sub   = "0.3"       # grey
    c_fill_pos = "#2ca02c"  # green (deposition)
    c_fill_neg = "#ff7f0e"  # orange (erosion)
    c_swl   = "#00aadd"   # cyan-blue (SWL)
    c_twl   = "#aa00dd"   # purple (TWL)

    # Raw dict insertion-order arrays – used for the profile line plots so that
    # the saved PNG matches plt.plot(list(dune_x.values()), list(dune_z.values())).
    # The dict key order is always the correct physical landward→seaward ordering
    # of profile pivots and must never be changed by a sort.
    x_init = np.array(list(dune_x_init.values()))
    z_init = np.array(list(dune_z_init.values()))

    x_curr = np.array(list(dune_x_curr.values()))
    z_curr = np.array(list(dune_z_curr.values()))

    # Sorted arrays – used only for np.interp, which requires monotonically
    # increasing x.  Pivots such as 'crest_temp' and 'seaward_crest' can end up
    # at very close but non-equal x-values (e.g. 183.44 vs 183) after growth,
    # making the insertion order non-monotonic.  Sorting fixes that for
    # interpolation without touching the line-plot order above.
    order_i = np.argsort(x_init, kind='stable')  # 'stable' preserves insertion order for equal-x pivots
    order_c = np.argsort(x_curr, kind='stable')  # 'stable' preserves insertion order for equal-x pivots

    x_init_s = x_init[order_i]
    z_init_s = z_init[order_i]
    x_curr_s = x_curr[order_c]
    z_curr_s = z_curr[order_c]

    xg_start = max(x_init_s[1:].min(), x_curr_s[1:].min())
    xg_end   = min(x_init_s[1:].max(), x_curr_s[1:].max())
    if x_init_s[1:].size == 0 or x_curr_s[1:].size == 0 or xg_start >= xg_end:
        # Profiles do not overlap – skip this snapshot
        print(f"Profile snapshot at {progress_pct:.0f}% skipped (profiles do not overlap).")
        return None

    xg = np.linspace(xg_start, xg_end, 500)

    z_init_g = np.interp(xg, x_init_s, z_init_s)
    z_curr_g = np.interp(xg, x_curr_s, z_curr_s)

    fig, ax = plt.subplots(figsize=(10, 4.8))

    ax.fill_between(xg, z_init_g, z_curr_g,
                    where=(z_curr_g > z_init_g),
                    color=c_fill_pos,
                    interpolate=True, alpha=0.35, zorder=1, label="Deposition")

    ax.fill_between(xg, z_init_g, z_curr_g,
                    where=(z_curr_g < z_init_g),
                    color=c_fill_neg,
                    interpolate=True, alpha=0.35, zorder=1, label="Erosion")

    ax.plot(x_init[1:], z_init[1:], color=c_init, lw=2.8, zorder=3, label="Initial profile")
    ax.plot(x_curr[1:], z_curr[1:], color=c_final, lw=2.8, zorder=3,
            label=f"Profile at {progress_pct:.0f}%")
    # ---- Boundary conditions ----
    if bc_info is not None:
        swl = bc_info.get('swl')
        twl = bc_info.get('total_water_level')
        hs  = bc_info.get('wave_height')
        tp  = bc_info.get('wave_period')
        r   = bc_info.get('runup')
        wind_speed  = bc_info.get('wind_speed')
        wind_dir    = bc_info.get('wind_direction')

        x_ramp_toe  = dune_x_curr.get('ramp_toe', x_curr[1:].max())
        x_shoreline = dune_x_curr.get('shoreline', x_curr[1:].max())
        x_sea_end   = x_shoreline + 0.1 * (x_shoreline - x_ramp_toe)

        # Extended profile grid to include the seaward visual margin
        _xfill = np.concatenate([xg, [x_sea_end]])
        _zfill = np.concatenate([z_curr_g, [0.0]])

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

        if swl is not None:
            # Find where the profile touches SWL (seaward-most point) then fill
            # only the seaward water body (no fill on the lee side of the dike).
            _tx_swl = _find_touch_x(_xfill, _zfill, swl)
            ax.fill_between(_xfill, swl, _zfill,
                            where=(_zfill <= swl) & (_xfill >= _tx_swl),
                            color=c_swl, alpha=0.25, interpolate=True,
                            zorder=1, label=f"SWL = {swl:.2f} m")
            ax.plot([_tx_swl, x_sea_end], [swl, swl],
                    color=c_swl, lw=2.0, linestyle='-.', zorder=5)

        if twl is not None:
            _dune_crest = dune_z_curr.get('seaward_crest')
            _is_overwash = (_dune_crest is not None and twl >= _dune_crest)
            if _is_overwash:
                # Overwash: fill from x=0 all the way seaward
                _z_back = dune_z_curr.get('back_ref',
                                          z_curr_g[0] if len(z_curr_g) > 0 else 0.0)
                _xfill_ow = np.concatenate([[0.0], xg, [x_sea_end]])
                _zfill_ow = np.concatenate([[_z_back], z_curr_g, [0.0]])
                ax.fill_between(_xfill_ow, twl, _zfill_ow,
                                where=(_zfill_ow <= twl),
                                color=c_twl, alpha=0.20, interpolate=True,
                                zorder=1, label=f"TWL = {twl:.2f} m (overwash)")
                ax.plot([0.0, x_sea_end], [twl, twl],
                        color=c_twl, lw=1.8, linestyle=':', zorder=5)
            else:
                # Find touch point then restrict fill to the seaward side only.
                _tx_twl = _find_touch_x(_xfill, _zfill, twl)
                ax.fill_between(_xfill, twl, _zfill,
                                where=(_zfill <= twl) & (_xfill >= _tx_twl),
                                color=c_twl, alpha=0.20, interpolate=True,
                                zorder=1, label=f"TWL = {twl:.2f} m")
                ax.plot([_tx_twl, x_sea_end], [twl, twl],
                        color=c_twl, lw=1.8, linestyle=':', zorder=5)

        # BC annotation box (upper-right)
        bc_lines = []
        if hs         is not None: bc_lines.append(f"Hs   = {hs:.2f} m")
        if tp         is not None: bc_lines.append(f"Tp   = {tp:.1f} s")
        if wind_speed is not None: bc_lines.append(f"U10  = {wind_speed:.1f} m/s")
        if wind_dir   is not None: bc_lines.append(f"Wdir = {wind_dir:.0f}\u00b0")
        if swl        is not None: bc_lines.append(f"SWL  = {swl:.2f} m")
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
    ax.legend(loc="lower right", fontsize=11, frameon=True)

    fig.tight_layout()

    plot_filename = f"{output_base}_{timestep}.png"
    fig.savefig(plot_filename, dpi=150)
    plt.close(fig)

    print(f"Profile snapshot saved to {plot_filename}")
    return plot_filename


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


def run_model(config_path, check_mass_conservation=False,
              plot_evolution=False, plot_evolution_interval=5.0,
              save_evolution_movie=False):
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

    volat = 0
    volab = 0
    volad = 0
    volrb = 0

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
            if r + swl <= dune_z['seaward_crest']:
                qd = vegred * 4 * p.cimpact * (r + swl - p.dfoot) ** 2 / tp
            else:
                print('OVERWASH')
                qd = vegred * 4 * p.cimpact * (((r + swl - p.dfoot) * (dune_z['seaward_crest']-p.dfoot)) / tp)
                alfa = (((r + swl - p.dfoot) / (dune_z['seaward_crest']-p.dfoot)) - 1) / aover
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
        qwnet = qwnet * max(1.0 - 4.2 * (dune_z['seaward_crest']-p.dfoot) / (dune_x['shoreline'] - dune_x['ramp_toe']), 0.0)

        ######## AEOLIAN DEPOSITION AT DUNE ########

        #if i>73596:
            #if dune_x['crest_temp']>=dune_x['seaward_crest']:
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
        vdune = partial_volume(dune_x, dune_z, substrate_x, substrate_z, dune_x['lee_crest'], dune_x['ramp_toe'])
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
            'eroded_volume': deltavramp_ero
        }
        
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

            # ---- Erosion-event snapshots ----
            _is_erosive = (deltavramp_ero > 0)
            if _is_erosive and not _erosion_active:
                # First timestep of a new erosion event
                _erosion_active = True
                _erosion_event_steps = 1
                _erosion_mid_saved = False
                _fp = _save_profile_snapshot(
                    results['dune_x'][0], results['dune_z'][0],
                    dune_x, dune_z,
                    substrate_x, substrate_z,
                    _plot_output_base, i + 1, current_pct,
                    bc_info=_bc_info
                )
                if _fp is not None:
                    _plot_saved_frames.append(_fp)
            elif _is_erosive and _erosion_active:
                _erosion_event_steps += 1
                if not _erosion_mid_saved and _erosion_event_steps >= 10:
                    # Mid-event snapshot (after at least 10 erosive steps)
                    _erosion_mid_saved = True
                    _fp = _save_profile_snapshot(
                        results['dune_x'][0], results['dune_z'][0],
                        dune_x, dune_z,
                        substrate_x, substrate_z,
                        _plot_output_base, i + 1, current_pct,
                        bc_info=_bc_info
                    )
                    if _fp is not None:
                        _plot_saved_frames.append(_fp)
            elif not _is_erosive and _erosion_active:
                # Last timestep of the erosion event – save end-of-event snapshot
                _erosion_active = False
                _fp = _save_profile_snapshot(
                    results['dune_x'][0], results['dune_z'][0],
                    dune_x, dune_z,
                    substrate_x, substrate_z,
                    _plot_output_base, i + 1, current_pct,
                    bc_info=_bc_info
                )
                if _fp is not None:
                    _plot_saved_frames.append(_fp)

            # ---- Regular interval snapshots ----
            if current_pct >= _plot_next_pct:
                _frame_path = _save_profile_snapshot(
                    results['dune_x'][0], results['dune_z'][0],
                    dune_x, dune_z,
                    substrate_x, substrate_z,
                    _plot_output_base, i + 1, current_pct,
                    bc_info=_bc_info
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

    if plot_evolution and save_evolution_movie and _plot_saved_frames:
        movie_path = f"{_plot_output_base}_evolution.gif"
        _create_evolution_movie(_plot_saved_frames, movie_path)
    
    ############################################################################
    # MASS CONSERVATION CHECK (if requested)
    ############################################################################
    
    if check_mass_conservation:
        print()  # Empty line for better formatting
        from cs_model.mass_conservation import check_mass_conservation as check_mass
        mass_check_result = check_mass(results, config, verbose=True)
        
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
            print(f"Total Volume Error: {mass_check_result['total_error']:+.6f} m³/m ({mass_check_result['total_error_percent']:.6f}%)")
            print(f"Initial Total Volume: {mass_check_result['initial_total']:.4f} m³/m")
            print(f"Final Total Volume: {mass_check_result['final_total']:.4f} m³/m")
            print(f"Error Threshold: {mass_check_result['error_threshold']:.1f}%")
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
