#!/usr/bin/env python
"""
Test that sediment conservation is correctly calculated when a dike limits
hydrodynamic (wave) erosion.

When the dike blocks the full potential erosion, the reported eroded volume
(deltavramp_ero) must equal the actual geometric dune volume decrease. A bug
existed where ramp_dike_temp was not updated in the "erosion beyond lee_crest"
case, causing a non-monotonic dune profile and a large spurious conservation
error at the dune.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from cs_model.main import run_model
from cs_model.mass_conservation import calc_dune_vol


CONFIG = os.path.join(os.path.dirname(__file__), 'test_dike_conservation.yaml')


def test_dike_erosion_conservation():
    """
    Run a short simulation with the dune seaward crest at the dike seaward
    toe (ysprim == ys_hs) and a high still-water level (swl=10 m) that
    forces severe wave erosion blocked by the dike.

    Conservation requirement: at every timestep the actual dune-volume
    change must match the change predicted from the reported eroded volume
    (deltavramp_ero) and aeolian deposition fluxes.  The per-timestep dune
    conservation error must remain within a tight numerical tolerance.
    """
    results = run_model(CONFIG, check_mass_conservation=False)

    substrate_x = results['substrate_x']
    substrate_z = results['substrate_z']
    n = len(results['time'])

    max_dune_error = 0.0

    for i in range(n):
        # Volumes at this step and the previous step
        vol_cur = calc_dune_vol(
            results['dune_x'][i + 1], results['dune_z'][i + 1],
            substrate_x, substrate_z
        )
        vol_prev = calc_dune_vol(
            results['dune_x'][i], results['dune_z'][i],
            substrate_x, substrate_z
        )
        actual_change = vol_cur - vol_prev

        # Expected dune change: aeolian deposition minus wave erosion
        aeolian = (results['qs_ramp'][i]
                   + results['qs_crest'][i]
                   + results['qs_lee'][i])
        eroded = results['eroded_volume'][i]
        expected_change = aeolian - eroded

        error = abs(actual_change - expected_change)
        if error > max_dune_error:
            max_dune_error = error

    # Allow a small numerical tolerance (1e-4 m³/m per timestep)
    tolerance = 1e-4
    assert max_dune_error < tolerance, (
        f"Dune conservation error too large when dike blocks erosion: "
        f"max per-timestep error = {max_dune_error:.6e} m³/m "
        f"(tolerance {tolerance:.0e} m³/m). "
        f"ramp_dike_temp may not be updated correctly after "
        f"erosion beyond lee_crest."
    )


if __name__ == "__main__":
    try:
        test_dike_erosion_conservation()
        print("✓ test_dike_erosion_conservation PASSED")
        sys.exit(0)
    except AssertionError as e:
        print(f"✗ test_dike_erosion_conservation FAILED: {e}")
        sys.exit(1)
