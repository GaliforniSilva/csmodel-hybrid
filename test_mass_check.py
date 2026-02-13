"""
Script to check mass conservation in the CS model.
Tracks all sediment fluxes and volumes throughout the simulation.
"""
import sys
import os
import numpy as np

# Add src directory to path
sys.path.insert(0, '/home/runner/work/csmodel-hybrid/csmodel-hybrid/src')

from cs_model.main import run_model
from cs_model.functions import partial_volume, total_volume
import pickle

def check_mass_conservation(config_path):
    """
    Run the model and check mass conservation at each time step.
    """
    print("Running model with mass conservation checks...")
    print("=" * 80)
    
    # Run the model
    run_model(config_path)
    
    # Load results
    output_filename = 'cs_model_test6.pkl'
    with open(output_filename, 'rb') as f:
        results = pickle.load(f)
    
    print("\nAnalyzing mass conservation...")
    print("=" * 80)
    
    # Calculate initial and final volumes
    initial_dune_x = results['dune_x'][0]
    initial_dune_z = results['dune_z'][0]
    substrate_x = results['substrate_x']
    substrate_z = results['substrate_z']
    
    final_dune_x = results['dune_x'][-1]
    final_dune_z = results['dune_z'][-1]
    
    # Calculate total volumes
    initial_total = total_volume(
        [initial_dune_x[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                      'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                      'ramp_toe_temp', 'ramp_toe', 'shoreline']],
        [initial_dune_z[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                      'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                      'ramp_toe_temp', 'ramp_toe', 'shoreline']],
        substrate_x, substrate_z
    )
    
    final_total = total_volume(
        [final_dune_x[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                    'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                    'ramp_toe_temp', 'ramp_toe', 'shoreline']],
        [final_dune_z[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                    'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                    'ramp_toe_temp', 'ramp_toe', 'shoreline']],
        substrate_x, substrate_z
    )
    
    # Sum all fluxes
    total_qs_ramp = np.sum(results['qs_ramp'])
    total_qs_crest = np.sum(results['qs_crest'])
    total_qs_lee = np.sum(results['qs_lee'])
    total_qwnet = np.sum(results['qwnet'])
    total_qbar = np.sum(results['qbar'])
    total_eroded = np.sum(results['eroded_volume'])
    total_overwash = np.sum(results['overwash_volume'])
    
    # Calculate volume changes
    initial_vbeach = results['vbeach'][0]
    final_vbeach = results['vbeach'][-1]
    initial_volbar = results['volbar'][0]
    final_volbar = results['volbar'][-1]
    
    print(f"\nINITIAL STATE:")
    print(f"  Total dune volume: {initial_total:.4f} m³/m")
    print(f"  Beach volume: {initial_vbeach:.4f} m³/m")
    print(f"  Bar volume: {initial_volbar:.4f} m³/m")
    print(f"  Total initial: {initial_total + initial_vbeach + initial_volbar:.4f} m³/m")
    
    print(f"\nFINAL STATE:")
    print(f"  Total dune volume: {final_total:.4f} m³/m")
    print(f"  Beach volume: {final_vbeach:.4f} m³/m")
    print(f"  Bar volume: {final_volbar:.4f} m³/m")
    print(f"  Total final: {final_total + final_vbeach + final_volbar:.4f} m³/m")
    
    print(f"\nVOLUME CHANGES:")
    print(f"  Dune volume change: {final_total - initial_total:.4f} m³/m")
    print(f"  Beach volume change: {final_vbeach - initial_vbeach:.4f} m³/m")
    print(f"  Bar volume change: {final_volbar - initial_volbar:.4f} m³/m")
    print(f"  Net volume change: {(final_total + final_vbeach + final_volbar) - (initial_total + initial_vbeach + initial_volbar):.4f} m³/m")
    
    print(f"\nFLUXES (integrated over time):")
    print(f"  Aeolian deposition on ramp: {total_qs_ramp:.4f} m³/m")
    print(f"  Aeolian deposition on crest: {total_qs_crest:.4f} m³/m")
    print(f"  Aeolian deposition on lee: {total_qs_lee:.4f} m³/m")
    print(f"  Net aeolian flux from beach: {-total_qwnet:.4f} m³/m")
    print(f"  Bar-beach exchange: {total_qbar:.4f} m³/m")
    print(f"  Total eroded from dune: {-total_eroded:.4f} m³/m")
    print(f"  Total overwash volume: {total_overwash:.4f} m³/m")
    
    # Mass balance check
    print(f"\nMASS BALANCE CHECK:")
    expected_dune_change = total_qs_ramp + total_qs_crest + total_qs_lee - total_eroded
    actual_dune_change = final_total - initial_total
    dune_balance_error = actual_dune_change - expected_dune_change
    
    print(f"  Expected dune change (from fluxes): {expected_dune_change:.4f} m³/m")
    print(f"  Actual dune change: {actual_dune_change:.4f} m³/m")
    print(f"  Dune balance error: {dune_balance_error:.4f} m³/m")
    print(f"  Dune balance error (%): {100 * dune_balance_error / (initial_total + 1e-10):.2f}%")
    
    expected_beach_change = -total_qwnet + total_qbar - (total_qs_ramp + total_qs_crest + total_qs_lee) + total_eroded
    actual_beach_change = final_vbeach - initial_vbeach
    beach_balance_error = actual_beach_change - expected_beach_change
    
    print(f"\n  Expected beach change (from fluxes): {expected_beach_change:.4f} m³/m")
    print(f"  Actual beach change: {actual_beach_change:.4f} m³/m")
    print(f"  Beach balance error: {beach_balance_error:.4f} m³/m")
    print(f"  Beach balance error (%): {100 * beach_balance_error / (initial_vbeach + 1e-10):.2f}%")
    
    total_balance_error = dune_balance_error + beach_balance_error
    print(f"\n  TOTAL MASS BALANCE ERROR: {total_balance_error:.4f} m³/m")
    print(f"  Total error (%): {100 * total_balance_error / (initial_total + initial_vbeach + 1e-10):.2f}%")
    
    if abs(total_balance_error) > 0.01:  # 1 cm³/m threshold
        print("\n⚠️  WARNING: Significant mass balance error detected!")
        return False
    else:
        print("\n✓ Mass is approximately conserved.")
        return True

if __name__ == "__main__":
    config_path = '/home/runner/work/csmodel-hybrid/csmodel-hybrid/test_mass_conservation.yaml'
    check_mass_conservation(config_path)
