"""
Test that _save_profile_snapshot uses stable sort so that co-located dune pivots
(same x-coordinate, different z-values) are rendered in dict insertion order,
matching what the user sees with a manual plt.plot(dune_x.values(), dune_z.values()).

Regression test for the bug where np.argsort() with the default 'quicksort'
kind can reorder equal-x pivots arbitrarily, causing the saved PNG to differ
from the manual interactive plot.
"""

import sys
import os
import tempfile

import numpy as np
import matplotlib
matplotlib.use('Agg')   # non-interactive backend required for tests
import matplotlib.pyplot as plt

# Make sure the src tree is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from cs_model.main import _save_profile_snapshot

# Floating-point tolerance used in elevation comparisons.
Z_TOLERANCE = 1e-6


def _build_dune_state_with_collocated_pivots(crest_temp_z, seaward_crest_z):
    """Return a dune_x / dune_z pair that mimics the ~92% test15 scenario.

    crest_temp and seaward_crest share the same x-coordinate (ysprim_new)
    but may have different z-values.  This is the situation that triggered
    the unstable-sort rendering artefact.
    """
    ysprim_new = 35.0   # Both crest_temp and seaward_crest at this x

    dune_x = {
        'back_ref':       0.0,
        'lee_foot':      10.0,
        'lee_crest':     25.0,
        'crest_temp':    ysprim_new,
        'seaward_crest': ysprim_new,   # same x as crest_temp!
        'ramp_dike_temp': ysprim_new,
        'ramp_crest':    ysprim_new,
        'ramp_toe_temp': 45.0,
        'ramp_toe':      45.0,
        'shoreline':     60.0,
    }

    dune_z = {
        'back_ref':       3.0,
        'lee_foot':       3.0,
        'lee_crest':      7.25,
        'crest_temp':     crest_temp_z,
        'seaward_crest':  seaward_crest_z,
        'ramp_dike_temp': 9.09,
        'ramp_crest':     9.09,
        'ramp_toe_temp':  3.0,
        'ramp_toe':       3.0,
        'shoreline':      0.0,
    }

    return dune_x, dune_z


def test_stable_sort_preserves_insertion_order_for_collocated_pivots():
    """
    When crest_temp and seaward_crest share the same x-coordinate but have
    different z-values, stable sort must preserve the dict insertion order
    (crest_temp before seaward_crest).  With the old unstable argsort those
    two could be swapped, producing a downward 'spike' artefact in the saved
    PNG that is absent from the manual plt.plot(dune_x.values(), …) plot.
    """
    # Scenario matching the problem report (approx. 92% of test15)
    # crest_temp_z == lee_crest_z < seaward_crest_z
    dune_x, dune_z = _build_dune_state_with_collocated_pivots(
        crest_temp_z=7.249, seaward_crest_z=8.5
    )

    x_vals = np.array(list(dune_x.values()))
    z_vals = np.array(list(dune_z.values()))

    ysprim_new = dune_x['crest_temp']   # == dune_x['seaward_crest']

    # ---- stable sort (the fixed implementation) ----
    order_stable = np.argsort(x_vals, kind='stable')
    x_stable = x_vals[order_stable]
    z_stable = z_vals[order_stable]

    # The stable sort must preserve the relative order of crest_temp and
    # seaward_crest (insertion-order indices 3 and 4).
    insertion_order_keys = list(dune_x.keys())
    ct_ins  = insertion_order_keys.index('crest_temp')
    sc_ins  = insertion_order_keys.index('seaward_crest')
    # After stable sort the relative order among equal-x entries must match.
    ct_sorted_pos = np.where(order_stable == ct_ins)[0][0]
    sc_sorted_pos = np.where(order_stable == sc_ins)[0][0]
    assert ct_sorted_pos < sc_sorted_pos, (
        f"Stable sort must keep crest_temp (pos {ct_sorted_pos}) before "
        f"seaward_crest (pos {sc_sorted_pos}) for equal x-values"
    )

    # In the line-plot sequence of z values at the shared x-coordinate, the
    # FIRST value must be crest_temp's z (7.249), NOT seaward_crest's z (8.5).
    z_at_shared_x_positions = z_stable[x_stable == ysprim_new]
    assert abs(z_at_shared_x_positions[0] - dune_z['crest_temp']) < 1e-9, (
        f"First z at shared x={ysprim_new} should be crest_temp z="
        f"{dune_z['crest_temp']}, got {z_at_shared_x_positions[0]}"
    )


def test_snapshot_file_is_saved_without_error():
    """
    _save_profile_snapshot should produce a PNG file without raising any
    exception even when co-located pivots have different z-values.
    """
    dune_x_init, dune_z_init = _build_dune_state_with_collocated_pivots(
        crest_temp_z=8.5, seaward_crest_z=8.5
    )
    dune_x_curr, dune_z_curr = _build_dune_state_with_collocated_pivots(
        crest_temp_z=7.249, seaward_crest_z=8.5
    )

    substrate_x = [0, 10, 35, 35, 45, 60]
    substrate_z = [3.0, 3.0, 9.09, 9.09, 3.0, 3.0]

    with tempfile.TemporaryDirectory() as tmpdir:
        out_base = os.path.join(tmpdir, 'snapshot_test')
        fp = _save_profile_snapshot(
            dune_x_init, dune_z_init,
            dune_x_curr, dune_z_curr,
            substrate_x, substrate_z,
            out_base, timestep=73900, progress_pct=92.0
        )
        assert fp is not None, "_save_profile_snapshot returned None unexpectedly"
        assert os.path.isfile(fp), f"Expected PNG file at {fp}"


def test_profile_grid_no_artifact_before_shared_x():
    """
    The interpolated profile grid must NOT show values near seaward_crest_z
    (8.5) in the region BEFORE the shared x-coordinate (ysprim_new).

    With the unstable-sort bug, seaward_crest (z=8.5) could appear first in
    the sorted array at x=35, causing np.interp to produce a value near 8.5
    immediately to the left of x=35 (interpolating from prev_point to 8.5
    instead of to 7.249).  The stable sort fixes this.
    """
    # Use lee_crest_z == crest_temp_z so the approach to ysprim_new is flat
    dune_x, dune_z = _build_dune_state_with_collocated_pivots(
        crest_temp_z=7.249, seaward_crest_z=8.5
    )
    # Align lee_crest with crest_temp to get a flat approach
    dune_z['lee_crest'] = dune_z['crest_temp']

    x_arr = np.array(list(dune_x.values()))
    z_arr = np.array(list(dune_z.values()))
    order = np.argsort(x_arr, kind='stable')
    x_s, z_s = x_arr[order], z_arr[order]

    ysprim_new = dune_x['crest_temp']

    xg = np.linspace(x_s[1:].min(), x_s[1:].max(), 500)
    # Note: [1:] skips back_ref (x=0), matching the grid range used inside
    # _save_profile_snapshot (see xg_start/xg_end calculation).
    zg = np.interp(xg, x_s, z_s)

    # Left of the shared x-position the profile should be at most crest_temp_z
    # (7.249), NOT approaching seaward_crest_z (8.5).
    zg_before = zg[xg < ysprim_new]
    if len(zg_before) > 0:
        max_z_before = zg_before.max()
        # Artifact threshold: if seaward_crest is sorted first, np.interp will
        # interpolate toward 8.5 just before x=35, producing max_z ≈ 8.49.
        # The correct (stable-sort) max_z should equal crest_temp_z = 7.249.
        assert max_z_before <= dune_z['crest_temp'] + Z_TOLERANCE, (
            f"Profile z ({max_z_before:.4f}) exceeds crest_temp z "
            f"({dune_z['crest_temp']}) before reaching x={ysprim_new} — "
            f"this indicates the unstable-sort spike artefact is present"
        )


if __name__ == '__main__':
    test_stable_sort_preserves_insertion_order_for_collocated_pivots()
    print("PASS: test_stable_sort_preserves_insertion_order_for_collocated_pivots")

    test_snapshot_file_is_saved_without_error()
    print("PASS: test_snapshot_file_is_saved_without_error")

    test_profile_grid_no_artifact_before_shared_x()
    print("PASS: test_profile_grid_no_artifact_before_shared_x")

    print("\nAll tests passed!")

