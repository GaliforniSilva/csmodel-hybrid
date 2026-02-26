#!/usr/bin/env python
"""
Test conservation of overwash volume when the dune crest migrates behind the dike.

Before the fix, the Delta Volume Overwash values would jump to ~+8.8 and then
to ~ -6.7 at the time step the dune fully crosses behind the dike, with
subsequent errors of -2 to -3.8 m²/m per step — a clear conservation violation.

The fix ensures that, when `seaward_crest <= ylprim_hs`, a volume-exact
numerical formulation is used for `f_crest_erosion` so that the sediment
removed from the crest equals exactly `erod_overwash`, which is then deposited
on the lee.  The seaward crest x-position is held fixed and `ramp_dike_temp`
is kept consistent with the new crest height on the dike back face.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from cs_model.main import run_model


def _collect_overwash_deltas(config_path):
    """Run the model and collect Delta Volume Overwash values from stdout."""
    import io
    from contextlib import redirect_stdout

    buf = io.StringIO()
    with redirect_stdout(buf):
        run_model(config_path)

    deltas = []
    for line in buf.getvalue().splitlines():
        if 'Delta Volume Overwash' in line:
            try:
                val = float(line.split('=')[-1].strip())
                deltas.append(val)
            except ValueError:
                pass
    return deltas


def test_overwash_conservation_dune_behind_dike():
    """
    Verify that |Delta Volume Overwash| is near machine precision for every
    time step, including and after the step when the dune crest crosses
    behind the dike.

    The test uses a high swl (10 m) so the dike is overtopped from the very
    first time step and the dune is pushed behind the dike early in the run.

    After the fix the lee update uses partial_volume in the brentq solver so
    the volume it converges to is exactly what partial_volume measures
    afterwards.  The threshold (1e-10 m²/m) reflects this machine-precision
    accuracy and detects any regression to the previous ~0.3 m²/m errors.
    """
    config = os.path.join(os.path.dirname(__file__), 'test_dune_behind_dike.yaml')
    deltas = _collect_overwash_deltas(config)

    assert len(deltas) > 0, "No overwash steps found in the simulation output"

    max_abs_delta = max(abs(d) for d in deltas)
    assert max_abs_delta < 1e-10, (
        f"Overwash conservation violation detected: "
        f"max |Delta Volume Overwash| = {max_abs_delta:.2e} m²/m "
        f"(threshold = 1e-10 m²/m). "
        f"All deltas: {[f'{d:.2e}' for d in deltas]}"
    )
    print(f"✓ Max |Delta Volume Overwash| = {max_abs_delta:.2e} m²/m  (all {len(deltas)} steps)")


def test_overwash_elevated_back_substrate():
    """
    Regression test for crash when base_elevation_back > dfoot.

    When the substrate behind the dike is higher than the dune foot elevation
    (base_elevation_back=6.0 > dfoot=3.0), the substrate intersection in the
    lee update raises the trial lee-foot elevation above lee_foot_z_cur.  This
    made the trial volume at the upper brentq bound exceed
    vdune_test_before_overwash, causing both bounds to share the same sign and
    brentq to raise ValueError('f(a) and f(b) must have different signs').

    The fix computes a consistent trial baseline at the upper bound so that
    brentq always receives a valid bracket.

    Note: the first overwash step may produce a non-negligible Delta Volume
    because the initial profile has the lee-foot below the substrate (an
    initialisation inconsistency).  Subsequent steps should be at machine
    precision.
    """
    config = os.path.join(os.path.dirname(__file__), 'test_katwijk_elevated_back.yaml')
    # The primary assertion is that no exception is raised (previously crashed)
    deltas = _collect_overwash_deltas(config)
    assert len(deltas) > 0, "No overwash steps found in the simulation output"

    # After the first step (which may have a large delta from the initial
    # substrate/foot inconsistency), all subsequent steps must conserve volume
    # at machine precision.
    if len(deltas) > 1:
        max_subsequent_delta = max(abs(d) for d in deltas[1:])
        assert max_subsequent_delta < 1e-10, (
            f"Overwash conservation violation after step 1: "
            f"max |Delta Volume Overwash| = {max_subsequent_delta:.2e} m²/m"
        )
    print(f"✓ Elevated-back-substrate case completed with {len(deltas)} overwash steps")
    print(f"  Step-1 delta = {deltas[0]:.2e} m²/m (initial substrate correction)")
    if len(deltas) > 1:
        print(f"  Max |delta| steps 2+: {max(abs(d) for d in deltas[1:]):.2e} m²/m")


if __name__ == '__main__':
    test_overwash_conservation_dune_behind_dike()
    test_overwash_elevated_back_substrate()
    print("All tests passed.")
