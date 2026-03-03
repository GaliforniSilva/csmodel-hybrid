#!/usr/bin/env python
"""
Test that ramp toe position is computed consistently with ramp_volume after
wave erosion moves ramp_toe inward without changing ramp_crest.

Before the fix, `morpho_update` used:
    ys_new_r = dune_x['ramp_crest'] + (sramp_new_r / betas)
which is only consistent with `ramp_volume` when the ramp geometry is
perfectly triangular (ramp_toe == ramp_crest + sramp/betas). After partial
wave erosion this condition no longer holds and the formula over-estimates
the new ramp toe position, depositing more volume than `qs_ramp` onto the
dune.

The fix restores the commented-out formula:
    ys_new_r = ys + ((sramp_new_r - sramp) / betas)
which matches what `ramp_volume` uses internally (yspotential = ys + ...).
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from cs_model.functions import morpho_update, total_volume, ramp_volume


def _dune_vol(dune_x, dune_z, substrate_x, substrate_z):
    """Calculate total dune volume above substrate."""
    order = ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp',
             'seaward_crest', 'ramp_dike_temp', 'ramp_crest',
             'ramp_toe_temp', 'ramp_toe', 'shoreline']
    x_list = [dune_x[k] for k in order if k in dune_x and dune_x[k] is not None]
    z_list = [dune_z[k] for k in order if k in dune_z and dune_z[k] is not None]
    return total_volume(x_list, z_list, substrate_x, substrate_z)


def test_ramp_toe_conservation_after_erosion():
    """
    Directly test the ramp toe formula with a geometrically inconsistent
    ramp (ramp_toe not at ramp_crest + sramp/betas, as occurs after
    partial wave erosion).  Aeolian accretion on the ramp must deposit
    exactly qs_ramp m³/m - no more, no less.
    """
    dfoot = 3.0
    betas = 0.1
    betal = 0.1
    s = 7.0     # crest height above dfoot
    sramp = 6.0  # ramp height above dfoot (< s, so brentq bounds are valid)

    # Substrate: flat beach below dune (all at z=3 seaward of the dike)
    substrate_x = [0, 72, 83, 189, 200, 500]
    substrate_z = [3.0, 3.0, 8.5, 8.5, 3.0, 3.0]

    # Build a dune profile where ramp_toe has been moved INWARD by prior
    # erosion but ramp_crest is unchanged.  Specifically:
    #   Normal toe position: ramp_crest + sramp/betas = 260 + 60 = 320
    #   After partial erosion: ramp_toe moved to 290 (< 320)
    # This is geometrically inconsistent and exercises the bug in the
    # original ys_new_r = dune_x['ramp_crest'] + (sramp_new_r / betas).
    ramp_crest_x = 260.0
    ramp_toe_eroded_x = 290.0  # moved inward by 30 m from original 320

    dune_x = {
        'back_ref': 0.0,
        'lee_foot': 120.0,
        'lee_crest': 180.0,
        'crest_temp': ramp_crest_x,
        'seaward_crest': ramp_crest_x,
        'ramp_dike_temp': ramp_crest_x,
        'ramp_crest': ramp_crest_x,
        'ramp_toe_temp': ramp_toe_eroded_x,
        'ramp_toe': ramp_toe_eroded_x,
        'shoreline': 450.0,
    }
    dune_z = {
        'back_ref': dfoot,
        'lee_foot': dfoot,
        'lee_crest': s + dfoot,
        'crest_temp': s + dfoot,
        'seaward_crest': s + dfoot,   # height = s + dfoot = 10.0
        'ramp_dike_temp': sramp + dfoot,
        'ramp_crest': sramp + dfoot,  # height = sramp + dfoot = 9.0 (< seaward_crest)
        'ramp_toe_temp': dfoot,
        'ramp_toe': dfoot,
        'shoreline': 0.0,
    }

    # Verify brentq bounds are distinct (required for the test to exercise the fix)
    lower = dune_z['ramp_crest'] - dfoot
    upper = dune_z['seaward_crest'] - dfoot
    assert lower < upper, (
        f"Test setup error: ramp_crest height ({lower}) must be < "
        f"seaward_crest height ({upper}) for brentq to have valid bounds."
    )

    # Verify geometric inconsistency (ramp_toe is NOT at ramp_crest + sramp/betas)
    expected_toe = ramp_crest_x + sramp / betas  # = 260 + 60 = 320
    assert abs(ramp_toe_eroded_x - expected_toe) > 1.0, (
        "Test setup error: ramp geometry should be inconsistent "
        f"(ramp_toe={ramp_toe_eroded_x} should differ from expected={expected_toe})"
    )

    dt = 600.0
    qs_ramp_deposit = 2.0  # m³/m deposited on the ramp

    vol_before = _dune_vol(dune_x, dune_z, substrate_x, substrate_z)

    dune_x2, dune_z2, _, _, _, _ = morpho_update(
        dune_x.copy(), dune_z.copy(),
        substrate_x, substrate_z,
        qs_ramp=qs_ramp_deposit, qs_crest=0.0, qs_lee=0.0,
        qs=0.0, ql=0.0,
        r=0.0, swl=0.0,  # no wave erosion
        vnourd=0.0, dt=dt,
        betas=betas, betas_hs=0.5,
        betal=betal, betal_hs=0.5,
        vramp=0.0, deltavall=0.0,
    )

    vol_after = _dune_vol(dune_x2, dune_z2, substrate_x, substrate_z)

    dune_vol_change = vol_after - vol_before
    conservation_error = abs(dune_vol_change - qs_ramp_deposit)

    assert conservation_error < 1e-6, (
        f"Dune volume change ({dune_vol_change:.8f} m³/m) should equal "
        f"qs_ramp ({qs_ramp_deposit:.8f} m³/m). "
        f"Conservation error: {conservation_error:.2e} m³/m. "
        f"ramp_toe before={ramp_toe_eroded_x}, ramp_crest={ramp_crest_x}, "
        f"sramp/betas={sramp/betas:.1f} (geometric inconsistency present)."
    )

