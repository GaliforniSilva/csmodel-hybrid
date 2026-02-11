"""
Comprehensive mass conservation report for CS Model.
"""
import sys
import numpy as np
import pickle
import yaml

sys.path.insert(0, 'src')
from cs_model.main import run_model
from cs_model.functions import total_volume

def calculate_actual_beach_volume(dune_x, vbeach, dclos, dfoot):
    """Calculate the actual beach volume (vbeachtot) from geometry."""
    shoreline = dune_x['shoreline']
    ramp_toe = dune_x['ramp_toe']
    db = vbeach / (shoreline - ramp_toe)
    vbeachtot = ramp_toe * (dclos + dfoot) + (shoreline - ramp_toe) * (dclos + db)
    return vbeachtot

def calc_dune_vol(dune_x, dune_z, substrate_x, substrate_z):
    """Calculate total dune volume above substrate."""
    x_list = [dune_x[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                   'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                   'ramp_toe_temp', 'ramp_toe', 'shoreline'] 
              if k in dune_x and dune_x[k] is not None]
    z_list = [dune_z[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                   'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                   'ramp_toe_temp', 'ramp_toe', 'shoreline']
              if k in dune_z and dune_z[k] is not None]
    return total_volume(x_list, z_list, substrate_x, substrate_z)

def check_mass_conservation(config_path):
    """
    Run the model and generate a comprehensive mass conservation report.
    """
    print("=" * 80)
    print("CS MODEL MASS CONSERVATION REPORT")
    print("=" * 80)
    print()
    
    # Load config to get parameters
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    dclos = config['erosion']['dclos']
    dfoot = config['dune_parameters']['dfoot']
    output_file = config.get('output', {}).get('filename', 'cs_model_results.pkl')
    
    print(f"Configuration file: {config_path}")
    print(f"Output file: {output_file}")
    print()
    
    # Run the model
    print("Running model...")
    run_model(config_path)
    print("Model run completed.")
    print()
    
    # Load results
    with open(output_file, 'rb') as f:
        results = pickle.load(f)
    
    substrate_x = results['substrate_x']
    substrate_z = results['substrate_z']
    n_steps = len(results['time'])
    dt = results['time'][1] - results['time'][0] if n_steps > 1 else 0
    total_time_hours = results['time'][-1] / 3600 if n_steps > 0 else 0
    
    print(f"Simulation details:")
    print(f"  Number of time steps: {n_steps}")
    print(f"  Time step: {dt:.0f} seconds")
    print(f"  Total simulation time: {total_time_hours:.1f} hours ({total_time_hours/24:.1f} days)")
    print()
    
    # Calculate volumes at start and end
    initial_dune_vol = calc_dune_vol(results['dune_x'][0], results['dune_z'][0], substrate_x, substrate_z)
    initial_beach_vol = calculate_actual_beach_volume(results['dune_x'][0], results['vbeach'][0], dclos, dfoot)
    initial_bar_vol = results['volbar'][0]
    initial_total = initial_dune_vol + initial_beach_vol + initial_bar_vol
    
    final_dune_vol = calc_dune_vol(results['dune_x'][-1], results['dune_z'][-1], substrate_x, substrate_z)
    final_beach_vol = calculate_actual_beach_volume(results['dune_x'][-1], results['vbeach'][-1], dclos, dfoot)
    final_bar_vol = results['volbar'][-1]
    final_total = final_dune_vol + final_beach_vol + final_bar_vol
    
    print("INITIAL STATE:")
    print(f"  Dune volume:  {initial_dune_vol:>12.4f} m³/m")
    print(f"  Beach volume: {initial_beach_vol:>12.4f} m³/m")
    print(f"  Bar volume:   {initial_bar_vol:>12.4f} m³/m")
    print(f"  Total volume: {initial_total:>12.4f} m³/m")
    print()
    
    print("FINAL STATE:")
    print(f"  Dune volume:  {final_dune_vol:>12.4f} m³/m")
    print(f"  Beach volume: {final_beach_vol:>12.4f} m³/m")
    print(f"  Bar volume:   {final_bar_vol:>12.4f} m³/m")
    print(f"  Total volume: {final_total:>12.4f} m³/m")
    print()
    
    # Calculate changes
    delta_dune = final_dune_vol - initial_dune_vol
    delta_beach = final_beach_vol - initial_beach_vol
    delta_bar = final_bar_vol - initial_bar_vol
    delta_total = final_total - initial_total
    
    print("VOLUME CHANGES:")
    print(f"  Δ Dune:  {delta_dune:>+12.4f} m³/m ({100*delta_dune/initial_dune_vol:>+7.2f}%)")
    print(f"  Δ Beach: {delta_beach:>+12.4f} m³/m ({100*delta_beach/initial_beach_vol:>+7.2f}%)")
    print(f"  Δ Bar:   {delta_bar:>+12.4f} m³/m ({100*delta_bar/initial_bar_vol:>+7.2f}%)")
    print(f"  Δ Total: {delta_total:>+12.4f} m³/m ({100*delta_total/initial_total:>+7.2f}%)")
    print()
    
    # Calculate integrated fluxes
    total_qs_ramp = np.sum(results['qs_ramp'][:n_steps])
    total_qs_crest = np.sum(results['qs_crest'][:n_steps])
    total_qs_lee = np.sum(results['qs_lee'][:n_steps])
    total_aeolian = total_qs_ramp + total_qs_crest + total_qs_lee
    
    total_qwnet_vol = np.sum([results['qwnet'][i] * dt for i in range(n_steps)])
    total_qbar_vol = np.sum([results['qbar'][i] * dt for i in range(n_steps)])
    total_eroded = np.sum(results['eroded_volume'][:n_steps])
    
    print("INTEGRATED FLUXES (over entire simulation):")
    print(f"  Aeolian transport from beach to dune: {total_aeolian:>12.4f} m³/m")
    print(f"    - Deposited on ramp:  {total_qs_ramp:>12.4f} m³/m")
    print(f"    - Deposited on crest: {total_qs_crest:>12.4f} m³/m")
    print(f"    - Deposited on lee:   {total_qs_lee:>12.4f} m³/m")
    print(f"  Bar-beach exchange (to bar if +):     {total_qbar_vol:>12.4f} m³/m")
    print(f"  Wave erosion (from dune):             {-total_eroded:>12.4f} m³/m")
    print()
    
    # Mass balance assessment
    print("MASS BALANCE ASSESSMENT:")
    print()
    print("Expected changes (based on fluxes):")
    expected_dune_change = total_aeolian - total_eroded
    expected_beach_change = -total_aeolian - total_qbar_vol + total_eroded
    expected_bar_change = total_qbar_vol
    expected_total_change = 0.0  # No external sources/sinks
    
    print(f"  Expected Δ Dune:  {expected_dune_change:>+12.4f} m³/m")
    print(f"  Expected Δ Beach: {expected_beach_change:>+12.4f} m³/m")
    print(f"  Expected Δ Bar:   {expected_bar_change:>+12.4f} m³/m")
    print(f"  Expected Δ Total: {expected_total_change:>+12.4f} m³/m")
    print()
    
    print("Mass balance errors:")
    dune_error = delta_dune - expected_dune_change
    beach_error = delta_beach - expected_beach_change
    bar_error = delta_bar - expected_bar_change
    total_error = delta_total - expected_total_change
    
    print(f"  Dune error:  {dune_error:>+12.6f} m³/m ({100*abs(dune_error)/initial_dune_vol:.4f}%)")
    print(f"  Beach error: {beach_error:>+12.6f} m³/m ({100*abs(beach_error)/initial_beach_vol:.4f}%)")
    print(f"  Bar error:   {bar_error:>+12.6f} m³/m ({100*abs(bar_error)/initial_bar_vol:.4f}%)")
    print(f"  TOTAL ERROR: {total_error:>+12.6f} m³/m ({100*abs(total_error)/initial_total:.6f}%)")
    print()
    
    # Assessment
    error_threshold = 0.1  # 0.1% threshold
    if abs(100*total_error/initial_total) < error_threshold:
        print("✓ RESULT: Mass is conserved within acceptable limits!")
        print(f"  (Error < {error_threshold}% threshold)")
    else:
        print("⚠ WARNING: Mass balance error exceeds acceptable threshold!")
        print(f"  (Error > {error_threshold}% threshold)")
    
    print()
    print("=" * 80)
    
    return abs(100*total_error/initial_total) < error_threshold

if __name__ == "__main__":
    if len(sys.argv) < 2:
        config_path = 'test_mass_conservation.yaml'
    else:
        config_path = sys.argv[1]
    
    check_mass_conservation(config_path)
