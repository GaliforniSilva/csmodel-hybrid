#!/usr/bin/env python
"""
Test script to verify conservation time-series tracking is working correctly.
"""

import sys
import os
import pickle

# Add src to path
sys.path.insert(0, 'src')

from cs_model.main import run_model

def test_conservation_timeseries():
    """Test that conservation time-series tracking works correctly."""
    
    print("=" * 80)
    print("Testing Conservation Time-Series Tracking")
    print("=" * 80)
    print()
    
    # Run a short test simulation
    print("Running test simulation...")
    results = run_model('test_short.yaml', check_mass_conservation=False)
    print()
    
    # Check that conservation data exists
    assert 'conservation' in results, "Conservation data not found in results!"
    print("✓ Conservation data found in results")
    
    conservation = results['conservation']
    
    # Check that all expected keys are present
    expected_keys = [
        'dune_error', 'beach_error', 'bar_error', 'total_error',
        'dune_error_pct', 'beach_error_pct', 'bar_error_pct', 'total_error_pct',
        'dune_volume', 'beach_volume', 'bar_volume', 'total_volume'
    ]
    
    for key in expected_keys:
        assert key in conservation, f"Missing key: {key}"
    print("✓ All expected conservation keys present")
    
    # Check that conservation data has correct length
    # Should be N+1 where N is number of timesteps (includes initial state)
    n_timesteps = len(results['time'])
    expected_length = n_timesteps + 1
    
    for key in expected_keys:
        actual_length = len(conservation[key])
        assert actual_length == expected_length, \
            f"Wrong length for {key}: expected {expected_length}, got {actual_length}"
    print(f"✓ Conservation data has correct length ({expected_length} points)")
    
    # Check that initial errors are zero
    assert conservation['dune_error'][0] == 0.0, "Initial dune error should be 0"
    assert conservation['beach_error'][0] == 0.0, "Initial beach error should be 0"
    assert conservation['bar_error'][0] == 0.0, "Initial bar error should be 0"
    assert conservation['total_error'][0] == 0.0, "Initial total error should be 0"
    print("✓ Initial errors are zero")
    
    # Check that initial volumes are positive
    assert conservation['dune_volume'][0] > 0, "Initial dune volume should be positive"
    assert conservation['beach_volume'][0] > 0, "Initial beach volume should be positive"
    assert conservation['bar_volume'][0] > 0, "Initial bar volume should be positive"
    assert conservation['total_volume'][0] > 0, "Initial total volume should be positive"
    print("✓ Initial volumes are positive")
    
    # Check that volumes change over time (at least one compartment)
    initial_total = conservation['total_volume'][0]
    final_total = conservation['total_volume'][-1]
    volume_changed = False
    for key in ['dune_volume', 'beach_volume', 'bar_volume']:
        if abs(conservation[key][-1] - conservation[key][0]) > 1e-6:
            volume_changed = True
            break
    assert volume_changed, "At least one compartment volume should change"
    print("✓ Volumes change over time")
    
    # Check that total volume is approximately conserved
    # (total error should be small relative to total volume)
    final_error = abs(conservation['total_error'][-1])
    relative_error = final_error / initial_total
    assert relative_error < 0.01, f"Total volume not conserved: {relative_error*100:.4f}% error"
    print(f"✓ Total volume approximately conserved (error: {relative_error*100:.6f}%)")
    
    # Print summary statistics
    print()
    print("Summary Statistics:")
    print(f"  Simulation timesteps: {n_timesteps}")
    print(f"  Conservation data points: {expected_length}")
    print(f"  Initial total volume: {initial_total:.4f} m³/m")
    print(f"  Final total volume: {final_total:.4f} m³/m")
    print(f"  Total volume change: {final_total - initial_total:.6f} m³/m")
    print(f"  Maximum dune error: {max(abs(e) for e in conservation['dune_error']):.6e} m³/m")
    print(f"  Maximum beach error: {max(abs(e) for e in conservation['beach_error']):.6e} m³/m")
    print(f"  Maximum bar error: {max(abs(e) for e in conservation['bar_error']):.6e} m³/m")
    print(f"  Maximum total error: {max(abs(e) for e in conservation['total_error']):.6e} m³/m")
    
    print()
    print("=" * 80)
    print("✓ All tests passed! Conservation time-series tracking is working correctly.")
    print("=" * 80)
    
    return True

if __name__ == "__main__":
    try:
        test_conservation_timeseries()
        sys.exit(0)
    except AssertionError as e:
        print()
        print("=" * 80)
        print("✗ Test failed!")
        print(f"Error: {e}")
        print("=" * 80)
        sys.exit(1)
    except Exception as e:
        print()
        print("=" * 80)
        print("✗ Unexpected error!")
        print(f"Error: {e}")
        print("=" * 80)
        import traceback
        traceback.print_exc()
        sys.exit(1)
