#!/usr/bin/env python
"""
Test that ramp_dike_temp z is consistent with seaward_crest z during recovery
when the dune crest is below the dike crest height and ramp_crest is on the
flat dike top (not at the exact ylprim_hs position).

Before the fix, line 1472 in functions.py set ramp_dike_temp z = s_hs + dfoot
(the full dike height = 8.5) regardless of seaward_crest z.  When
seaward_crest z < 8.5 this created a physically impossible vertical wall in
the profile at seaward_crest x, leading to the "misalignment in the detached
lee" reported in test15 overwash simulations.

The fix changes line 1472 to use dune_z['seaward_crest'] instead of
s_hs + dfoot, ensuring ramp_dike_temp is co-located (same z) with
seaward_crest in the scarp case.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

import numpy as np
from cs_model.functions import morpho_update


def test_ramp_dike_temp_consistent_with_seaward_crest_during_recovery():
    """
    Verify that ramp_dike_temp z equals seaward_crest z after a recovery step
    when:
      - seaward_crest is behind the dike (seaward_crest x < ylprim_hs)
      - seaward_crest z < dike crest height (s_hs + dfoot)
      - ramp_crest is on the flat dike top (ylprim_hs < ramp_crest < ysprim_hs)
      - qs_ramp > 0 (recovery / deposition, no erosion)

    In this scenario the recovery code reaches the 'else' (scarp) branch at
    what was line 1472.  Before the fix, ramp_dike_temp z was set to
    s_hs + dfoot = 8.5, creating a phantom vertical wall.  After the fix it
    equals seaward_crest z, keeping the profile geometrically consistent.
    """
    # Substrate: flat ground at z=3, dike from x=172 (z=3) to x=183 (z=8.5),
    # flat dike top 183–289 at z=8.5, seaward face to x=300 (z=3), beach to 500.
    substrate_x = [0.0, 172.0, 183.0, 289.0, 300.0, 500.0]
    substrate_z = [3.0,   3.0,   8.5,   8.5,   3.0,   3.0]

    # Dune state: ramp_crest is on the flat dike top (220, between 183 and 289)
    # seaward_crest is behind the dike (178 < 183) at height 7.2 < 8.5.
    # This represents the state at the end of an overwash period where the
    # erosion front stopped in the middle of the flat dike top.
    dune_x = {
        'back_ref':      0.0,
        'lee_foot':      110.0,
        'lee_crest':     150.0,
        'crest_temp':    178.0,
        'seaward_crest': 178.0,
        'ramp_dike_temp': 178.0,
        'ramp_crest':    220.0,   # on flat dike top (183 < 220 < 289)
        'ramp_toe_temp': 289.0,
        'ramp_toe':      300.0,
        'shoreline':     500.0,
    }
    dune_z = {
        'back_ref':      3.0,
        'lee_foot':      3.0,
        'lee_crest':     7.2,
        'crest_temp':    7.2,
        'seaward_crest': 7.2,     # below dike crest height (7.2 < 8.5)
        'ramp_dike_temp': 7.2,
        'ramp_crest':    8.5,     # on flat dike top at substrate height
        'ramp_toe_temp': 8.5,
        'ramp_toe':      3.0,
        'shoreline':     0.0,
    }

    # Parameters matching the test15 configuration
    betas    = 0.1
    betas_hs = 0.5
    betal    = 0.1
    betal_hs = 0.5
    vramp    = 0.0   # no existing ramp volume

    # Recovery step: small positive deposition on ramp, no erosion
    qs_ramp  = 5.0
    qs_crest = 0.0
    qs_lee   = 0.0
    qs       = 0.0   # no ramp erosion
    ql       = 0.0   # no overwash
    r        = 0.0   # no runup
    swl      = 0.0
    vnourd   = 0.0
    dt       = 600.0
    deltavall = 0.0

    morpho_update(
        dune_x, dune_z, substrate_x, substrate_z,
        qs_ramp, qs_crest, qs_lee, qs, ql,
        r, swl, vnourd, dt,
        betas, betas_hs, betal, betal_hs,
        vramp, deltavall,
    )

    # After the fix: ramp_dike_temp z must equal seaward_crest z.
    # Before the fix it was s_hs + dfoot = 8.5, creating a vertical wall.
    assert dune_z['ramp_dike_temp'] == dune_z['seaward_crest'], (
        f"ramp_dike_temp z ({dune_z['ramp_dike_temp']:.4f}) != "
        f"seaward_crest z ({dune_z['seaward_crest']:.4f}). "
        "A vertical wall was introduced at seaward_crest x, causing the "
        "detached-lee misalignment in overwash simulations."
    )

    # Additionally verify that ramp_dike_temp x == seaward_crest x
    # (they should be co-located in the scarp case).
    assert dune_x['ramp_dike_temp'] == dune_x['seaward_crest'], (
        f"ramp_dike_temp x ({dune_x['ramp_dike_temp']:.4f}) != "
        f"seaward_crest x ({dune_x['seaward_crest']:.4f})."
    )


if __name__ == '__main__':
    test_ramp_dike_temp_consistent_with_seaward_crest_during_recovery()
    print("All tests passed.")
