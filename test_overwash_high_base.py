#!/usr/bin/env python
"""
Test the dfoot_back parameter for cases where base_elevation_back > dfoot.

Without dfoot_back, dune_z['lee_foot'] is initialised to dfoot (3.0) while
the substrate at the back is base_elevation_back (6.0), placing the lee foot
below the substrate.  This creates an invalid brentq bracket during overwash
(ValueError: f(a) and f(b) must have different signs).

With dfoot_back: 6.0 in the config, dune_z['back_ref'] and dune_z['lee_foot']
are set to 6.0, making the geometry consistent and the model runs to completion.
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from cs_model.main import run_model


# Config with dfoot_back=6.0, base_elevation_back=6.0  (the corrected case)
_CONFIG_WITH_DFOOT_BACK = os.path.join(
    os.path.dirname(__file__), 'test_overwash_high_base.yaml'
)


def test_model_runs_with_dfoot_back():
    """
    Verify that setting dfoot_back >= base_elevation_back lets the model
    complete without a brentq ValueError.
    """
    results = run_model(_CONFIG_WITH_DFOOT_BACK)
    assert results is not None
    assert len(results['dune_x']) > 1


def test_dfoot_back_falls_back_when_absent():
    """
    Verify that omitting dfoot_back (existing configs) still works: the model
    initialises lee_foot to dfoot as before.
    """
    config = os.path.join(os.path.dirname(__file__), 'test_short.yaml')
    results = run_model(config)
    assert results is not None
    # lee_foot and back_ref should start at dfoot (3.0)
    assert results['dune_z'][0]['lee_foot'] == pytest.approx(3.0)
    assert results['dune_z'][0]['back_ref'] == pytest.approx(3.0)


def test_dfoot_back_sets_initial_lee_foot():
    """
    Verify that when dfoot_back is provided, dune_z['lee_foot'] and
    dune_z['back_ref'] are initialised to dfoot_back, not dfoot.
    """
    results = run_model(_CONFIG_WITH_DFOOT_BACK)
    # dfoot=3.0, dfoot_back=6.0 — initial lee_foot and back_ref must be 6.0
    assert results['dune_z'][0]['lee_foot'] == pytest.approx(6.0)
    assert results['dune_z'][0]['back_ref'] == pytest.approx(6.0)
    # ramp_toe should still use dfoot (3.0), not dfoot_back
    assert results['dune_z'][0]['ramp_toe'] == pytest.approx(3.0)
