"""
Test script to verify the ramp_crest height update fix.
Tests the f_ramp_dike_top case and checks mass conservation.
"""
import sys
import os
import numpy as np
import pickle

# Add src directory to path
sys.path.insert(0, '/home/runner/work/csmodel-hybrid/csmodel-hybrid/src')

from cs_model.main import run_model
from cs_model.functions import partial_volume, total_volume

def test_ramp_crest_fix():
    """
    Run the model with conditions that trigger f_ramp_dike_top and check:
    1. Dune growth rate is reasonable
    2. Mass is conserved
    """
    config_path = '/home/runner/work/csmodel-hybrid/csmodel-hybrid/test_ramp_crest_issue.yaml'
    
    print("=" * 80)
    print("TESTING RAMP_CREST HEIGHT UPDATE FIX")
    print("=" * 80)
    print(f"\nRunning model with config: {config_path}")
    print("This test uses parameters that trigger the f_ramp_dike_top case.")
    print()
    
    # Run the model
    run_model(config_path)
    
    # Load results
    output_filename = 'cs_model_test9b.pkl'
    with open(output_filename, 'rb') as f:
        results = pickle.load(f)
    
    print("\n" + "=" * 80)
    print("ANALYZING RESULTS")
    print("=" * 80)
    
    # Get key data
    substrate_x = results['substrate_x']
    substrate_z = results['substrate_z']
    times = results['time']
    
    # Calculate dune volumes over time
    dune_volumes = []
    ramp_crest_heights = []
    
    for i in range(len(results['dune_x'])):
        dune_x = results['dune_x'][i]
        dune_z = results['dune_z'][i]
        
        # Calculate total volume
        vol = total_volume(
            [dune_x[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                  'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                  'ramp_toe_temp', 'ramp_toe', 'shoreline']],
            [dune_z[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                  'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                  'ramp_toe_temp', 'ramp_toe', 'shoreline']],
            substrate_x, substrate_z
        )
        dune_volumes.append(vol)
        ramp_crest_heights.append(dune_z['ramp_crest'])
    
    # Print initial and final state
    print(f"\nINITIAL STATE:")
    print(f"  Dune volume: {dune_volumes[0]:.4f} m³/m")
    print(f"  Ramp crest height: {ramp_crest_heights[0]:.4f} m")
    print(f"  Ramp crest position: {results['dune_x'][0]['ramp_crest']:.4f} m")
    
    print(f"\nFINAL STATE (after {len(times)} time steps):")
    print(f"  Dune volume: {dune_volumes[-1]:.4f} m³/m")
    print(f"  Ramp crest height: {ramp_crest_heights[-1]:.4f} m")
    print(f"  Ramp crest position: {results['dune_x'][-1]['ramp_crest']:.4f} m")
    
    # Calculate changes
    volume_change = dune_volumes[-1] - dune_volumes[0]
    height_change = ramp_crest_heights[-1] - ramp_crest_heights[0]
    
    print(f"\nCHANGES:")
    print(f"  Dune volume change: {volume_change:.4f} m³/m")
    print(f"  Ramp crest height change: {height_change:.4f} m")
    
    # Sum all aeolian fluxes
    total_qs_ramp = np.sum(results['qs_ramp'])
    total_qs_crest = np.sum(results['qs_crest'])
    total_qs_lee = np.sum(results['qs_lee'])
    total_aeolian = total_qs_ramp + total_qs_crest + total_qs_lee
    
    # Sum erosion
    total_eroded = np.sum(results['eroded_volume'])
    
    print(f"\nSEDIMENT FLUXES:")
    print(f"  Total aeolian deposition: {total_aeolian:.4f} m³/m")
    print(f"  Total erosion: {total_eroded:.4f} m³/m")
    print(f"  Net flux to dune: {total_aeolian - total_eroded:.4f} m³/m")
    
    # Check mass conservation
    expected_change = total_aeolian - total_eroded
    mass_balance_error = volume_change - expected_change
    
    print(f"\nMASS CONSERVATION CHECK:")
    print(f"  Expected volume change (from fluxes): {expected_change:.4f} m³/m")
    print(f"  Actual volume change: {volume_change:.4f} m³/m")
    print(f"  Mass balance error: {mass_balance_error:.4f} m³/m")
    print(f"  Mass balance error (%): {100 * mass_balance_error / (dune_volumes[0] + 1e-10):.2f}%")
    
    # Check for reasonable growth
    print(f"\nGROWTH RATE CHECK:")
    if total_aeolian > 0:
        efficiency = volume_change / total_aeolian * 100
        print(f"  Deposition efficiency: {efficiency:.1f}%")
        print(f"  (Should be close to 100% if all aeolian sediment stays in dune)")
    
    # Determine if test passed
    print("\n" + "=" * 80)
    print("TEST RESULTS")
    print("=" * 80)
    
    passed = True
    
    # Check 1: Mass conservation (within 1%)
    if abs(mass_balance_error / (dune_volumes[0] + 1e-10)) > 0.01:
        print("❌ FAILED: Mass balance error exceeds 1%")
        passed = False
    else:
        print("✓ PASSED: Mass is conserved within 1%")
    
    # Check 2: Reasonable growth (dune should grow if there's net deposition)
    if expected_change > 0 and volume_change < 0.5 * expected_change:
        print("❌ FAILED: Dune growth is less than 50% of expected (indicates bug)")
        passed = False
    else:
        print("✓ PASSED: Dune growth rate is reasonable")
    
    # Check 3: Ramp crest height is reasonable (should be above dfoot)
    if ramp_crest_heights[-1] < substrate_z[4]:  # substrate_z[4] is dfoot
        print(f"❌ FAILED: Ramp crest height ({ramp_crest_heights[-1]:.2f} m) is below dfoot ({substrate_z[4]:.2f} m)")
        passed = False
    else:
        print(f"✓ PASSED: Ramp crest height is above dfoot")
    
    print("\n" + "=" * 80)
    if passed:
        print("✓ ALL TESTS PASSED")
    else:
        print("❌ SOME TESTS FAILED")
    print("=" * 80)
    
    return passed

if __name__ == "__main__":
    success = test_ramp_crest_fix()
    sys.exit(0 if success else 1)
