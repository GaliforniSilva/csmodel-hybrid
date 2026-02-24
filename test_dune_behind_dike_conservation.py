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
    Verify that |Delta Volume Overwash| is small for every time step,
    including and after the step when the dune crest crosses behind the dike.

    The test uses a high swl (10 m) so the dike is overtopped from the very
    first time step and the dune is pushed behind the dike early in the run.
    The threshold (0.5 m²/m) is generous enough to accept pre-existing minor
    numerical errors in other parts of the model while being tight enough to
    detect the ~8 m²/m jump that was present before the fix.
    """
    config = os.path.join(os.path.dirname(__file__), 'test_dune_behind_dike.yaml')
    deltas = _collect_overwash_deltas(config)

    assert len(deltas) > 0, "No overwash steps found in the simulation output"

    max_abs_delta = max(abs(d) for d in deltas)
    assert max_abs_delta < 0.5, (
        f"Large overwash conservation violation detected: "
        f"max |Delta Volume Overwash| = {max_abs_delta:.4f} m²/m "
        f"(threshold = 0.5 m²/m). "
        f"All deltas: {[round(d, 4) for d in deltas]}"
    )
    print(f"✓ Max |Delta Volume Overwash| = {max_abs_delta:.6f} m²/m  (all {len(deltas)} steps)")


if __name__ == '__main__':
    test_overwash_conservation_dune_behind_dike()
    print("All tests passed.")
