#!/usr/bin/env python
"""
Integration test for sediment conservation using the test15 example configuration.

The problem statement reported a large dune conservation error (~85 m³/m over
549 days, i.e. ~19.6% of initial dune volume) caused by an incorrect ramp toe
position formula in `morpho_update`.  This test runs the test15 configuration
for 800 hours (which is long enough to capture the first wave-erosion events
that trigger the bug) and asserts that mass is conserved within the model's
0.1% threshold.

Root cause fixed (functions.py):
    Before: ys_new_r = dune_x['ramp_crest'] + (sramp_new_r / betas)
    After:  ys_new_r = ys + ((sramp_new_r - sramp) / betas)
where `ys = dune_x['ramp_toe']` and `sramp = dune_z['ramp_crest'] - dfoot`.
The old formula was only correct when the ramp geometry was perfectly
triangular (ramp_toe == ramp_crest + sramp/betas).  After partial wave
erosion moves ramp_toe inward while leaving ramp_crest unchanged, the old
formula placed the new ramp toe too far seaward, depositing more sediment on
the dune than the aeolian flux warranted and creating a systematic positive
dune conservation error that accumulated over hundreds of erosion-accretion
cycles.
"""

import sys
import os
import tempfile

import yaml

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from cs_model.main import run_model
from cs_model.mass_conservation import check_mass_conservation

# Path to the test15 example configuration
TEST15_DIR = os.path.join(os.path.dirname(__file__),
                          'examples', 'sediment_conservation', 'test15')
TEST15_CONFIG = os.path.join(TEST15_DIR, 'config_test15.yaml')


def test_test15_sediment_conservation():
    """
    Run the test15 scenario for 800 hours and verify sediment conservation.

    800 hours is enough to include the first wave-erosion events (which begin
    around t = 744 h in the full run) that previously triggered the ramp-toe
    bug.  After the fix, the dune conservation error must be well within the
    accepted 0.1 % threshold.
    """
    with open(TEST15_CONFIG) as f:
        cfg = yaml.safe_load(f)

    # Shorten simulation to 800 h to keep CI runtime reasonable (~13 s)
    cfg['time_params']['time_length_hours'] = 800

    # Write the modified config to a temp file in the test15 directory so
    # that the relative paths to synthetic_swl.txt and synthetic_month.txt
    # resolve correctly.
    orig_dir = os.getcwd()
    os.chdir(TEST15_DIR)
    try:
        with tempfile.NamedTemporaryFile(
                mode='w', suffix='.yaml', dir='.', delete=False) as f:
            yaml.dump(cfg, f)
            tmp_cfg = f.name

        results = run_model(tmp_cfg, check_mass_conservation=False)
        mc = check_mass_conservation(results, cfg, verbose=False)
    finally:
        if os.path.exists(tmp_cfg):
            os.unlink(tmp_cfg)
        os.chdir(orig_dir)

    # The overall conservation check must pass (total error < 0.1 %)
    assert mc['conserved'], (
        f"test15 mass conservation FAILED after 800 h: "
        f"total error = {mc['total_error_percent']:.4f}% "
        f"(threshold = {mc['error_threshold']}%)"
    )

    # The dune specifically must have negligible error (< 0.1 m³/m absolute).
    # Before the fix this would have been ~5 m³/m for an 800 h run.
    assert abs(mc['dune_error']) < 0.1, (
        f"Dune conservation error too large after fix: "
        f"{mc['dune_error']:.4f} m³/m (should be < 0.1 m³/m). "
        f"The ramp-toe formula bug may have been reintroduced."
    )
