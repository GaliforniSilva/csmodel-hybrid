"""
test_flux_updates.py

Tests to verify that sediment fluxes are correctly updated when interacting
with the non-erodible layer (dike/substrate), ensuring proper sediment
conservation tracking.

Specifically tests:
1. erod_ramp is correctly set to 0 when non-erodible layer is exposed
   (deltavramp_ero = 0 returned from morpho_update)
2. erod_overwash is updated to the actual value (q_to_lee) when the
   non-erodible layer limits crest erosion during overwash
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from cs_model.functions import morpho_update


def make_base_profile():
    """Create a base dune/dike profile for testing."""
    # Dike: ylprim_hs=25, ysprim_hs=38, ys_hs=50, s_hs=2.0m
    substrate_x = [0.0, 10.0, 25.0, 38.0, 50.0, 100.0]
    substrate_z = [0.0, 0.0, 2.0, 2.0, 0.0, 0.0]  # dike height = 2.0m

    dfoot = 0.5  # dune foot elevation

    dune_x = {
        'back_ref': 0.0,
        'lee_foot': 15.0,
        'lee_crest': 22.0,
        'crest_temp': 30.0,
        'seaward_crest': 35.0,
        'ramp_dike_temp': 35.0,
        'ramp_crest': 35.0,
        'ramp_toe_temp': 45.0,
        'ramp_toe': 45.0,
        'shoreline': 60.0
    }

    dune_z = {
        'back_ref': dfoot,
        'lee_foot': dfoot,
        'lee_crest': 3.5,
        'crest_temp': 3.5,
        'seaward_crest': 3.5,
        'ramp_dike_temp': 2.0 + dfoot,
        'ramp_crest': 2.0 + dfoot,
        'ramp_toe_temp': dfoot,
        'ramp_toe': dfoot,
        'shoreline': 0.0
    }

    return dune_x, dune_z, substrate_x, substrate_z


def test_erod_ramp_zero_when_dike_exposed():
    """
    When the non-erodible layer is exposed, erod_ramp should be set to 0
    and deltavramp_ero should be 0 (no actual ramp erosion occurs).
    """
    dune_x, dune_z, substrate_x, substrate_z = make_base_profile()

    # Conditions that put runup on the dike (dike exposed, ramp can't erode)
    # r + swl is between dike foot (0.5) and dike crest (2.0+0.5=2.5)
    r = 1.5
    swl = 0.3
    qs = 0.01   # potential wave erosion flux > 0
    ql = 0.0    # no overwash
    dt = 60.0

    result = morpho_update(
        dune_x.copy(), dune_z.copy(), substrate_x, substrate_z,
        0.0, 0.0, 0.0, qs, ql, r, swl, 0.0, dt,
        0.25, 0.5, 0.25, 0.5, 0.0, 0.0
    )
    _, _, deltavramp, erod_overwash_out, deltavramp_ero, _ = result

    assert deltavramp_ero == 0.0, (
        f"Expected deltavramp_ero=0 when dike is exposed (no actual ramp erosion), "
        f"got {deltavramp_ero:.6f}"
    )
    assert erod_overwash_out == 0.0, f"Expected erod_overwash=0 (no overwash), got {erod_overwash_out}"
    print(f"  test_erod_ramp_zero_when_dike_exposed: PASS "
          f"(deltavramp_ero={deltavramp_ero:.6f}, erod_overwash={erod_overwash_out:.6f})")


def test_erod_overwash_updated_when_dike_limits():
    """
    When the non-erodible layer (dike) constrains crest erosion during overwash,
    erod_overwash returned from morpho_update must equal the ACTUAL volume
    transported to the lee side (q_to_lee), NOT the original potential value (ql*dt).

    This ensures results['overwash_volume'] records the actual overwash.
    """
    dune_x, dune_z, substrate_x, substrate_z = make_base_profile()

    # Large overwash: enough to lower crest below dike height (s_hs=2.0)
    # crest_width = seaward_crest - lee_crest = 35 - 22 = 13
    # V_cap(1.0m lowering) = 13*1 + 0.5*(1/0.25 + 1/0.25) = 13 + 4 = 17
    # erod_overwash = 30 > 17 -> dike will limit the erosion to ~17
    ql_potential = 30.0   # erod_overwash = ql * dt (dt=1) = 30
    dt = 1.0
    r = 3.2
    swl = 0.5  # r + swl = 3.7 > seaward_crest=3.5 (overwash) and > s_hs+dfoot=2.5 (dike submerged)

    result = morpho_update(
        dune_x.copy(), dune_z.copy(), substrate_x, substrate_z,
        0.0, 0.0, 0.0, 0.0, ql_potential, r, swl, 0.0, dt,
        0.25, 0.5, 0.25, 0.5, 0.0, 0.0
    )
    _, _, _, erod_overwash_out, deltavramp_ero, _ = result

    # The dike limits crest erosion: actual overwash < potential
    assert erod_overwash_out < ql_potential, (
        f"Expected erod_overwash ({erod_overwash_out:.4f}) < potential ({ql_potential:.4f}). "
        f"The non-erodible layer should limit overwash volume."
    )
    assert erod_overwash_out >= 0.0, "erod_overwash must be non-negative"
    print(f"  test_erod_overwash_updated_when_dike_limits: PASS "
          f"(erod_overwash={erod_overwash_out:.4f} < potential={ql_potential:.4f})")


def test_erod_overwash_unchanged_without_dike_limit():
    """
    When the dike does NOT limit crest erosion (small overwash), erod_overwash
    should remain equal to the potential value ql*dt.
    """
    dune_x, dune_z, substrate_x, substrate_z = make_base_profile()

    # Small overwash: V_cap(1.0) ≈ 17, use ql_potential=2 which won't trigger dike limit
    ql_potential = 2.0
    dt = 1.0
    r = 3.2
    swl = 0.5  # r + swl = 3.7 > seaward_crest=3.5 (overwash)

    result = morpho_update(
        dune_x.copy(), dune_z.copy(), substrate_x, substrate_z,
        0.0, 0.0, 0.0, 0.0, ql_potential, r, swl, 0.0, dt,
        0.25, 0.5, 0.25, 0.5, 0.0, 0.0
    )
    _, _, _, erod_overwash_out, _, _ = result

    # Without dike limit, erod_overwash should equal the potential value
    assert abs(erod_overwash_out - ql_potential) < 1e-10, (
        f"Expected erod_overwash ≈ {ql_potential:.4f} (no dike limit), "
        f"got {erod_overwash_out:.4f}"
    )
    print(f"  test_erod_overwash_unchanged_without_dike_limit: PASS "
          f"(erod_overwash={erod_overwash_out:.4f} == potential={ql_potential:.4f})")


def test_deltavramp_ero_correct_when_no_erosion():
    """
    When the non-erodible layer is fully exposed (erod_ramp=0),
    deltavramp_ero must be 0 - the flux going from dune to beach is 0.
    """
    dune_x, dune_z, substrate_x, substrate_z = make_base_profile()

    # Accretion case (no wave erosion): deltavramp_ero should be 0
    qs = 0.0
    ql = 0.0
    r = 0.3  # runup below dune foot
    swl = 0.1
    dt = 60.0

    result = morpho_update(
        dune_x.copy(), dune_z.copy(), substrate_x, substrate_z,
        0.005, 0.001, 0.0, qs, ql, r, swl, 0.0, dt,
        0.25, 0.5, 0.25, 0.5, 0.0, 0.0
    )
    _, _, _, erod_overwash_out, deltavramp_ero, _ = result

    assert deltavramp_ero == 0.0, (
        f"Expected deltavramp_ero=0 in accretion case, got {deltavramp_ero}"
    )
    assert erod_overwash_out == 0.0, (
        f"Expected erod_overwash=0 in accretion case, got {erod_overwash_out}"
    )
    print(f"  test_deltavramp_ero_correct_when_no_erosion: PASS "
          f"(deltavramp_ero={deltavramp_ero:.6f}, erod_overwash={erod_overwash_out:.6f})")


if __name__ == "__main__":
    print("=" * 70)
    print("Testing flux updates during dike/non-erodible layer interaction")
    print("=" * 70)
    print()

    passed = 0
    failed = 0
    tests = [
        test_erod_ramp_zero_when_dike_exposed,
        test_erod_overwash_updated_when_dike_limits,
        test_erod_overwash_unchanged_without_dike_limit,
        test_deltavramp_ero_correct_when_no_erosion,
    ]

    for test_fn in tests:
        try:
            test_fn()
            passed += 1
        except AssertionError as e:
            print(f"  {test_fn.__name__}: FAIL - {e}")
            failed += 1
        except Exception as e:
            print(f"  {test_fn.__name__}: ERROR - {e}")
            failed += 1

    print()
    print("=" * 70)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 70)

    sys.exit(0 if failed == 0 else 1)
