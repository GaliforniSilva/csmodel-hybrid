"""
Detailed mass conservation check for CS Model.
"""
import sys
import numpy as np
import pickle

sys.path.insert(0, 'src')
from cs_model.functions import total_volume, partial_volume

# Load results
with open('test_short.pkl', 'rb') as f:
    results = pickle.load(f)

substrate_x = results['substrate_x']
substrate_z = results['substrate_z']

# Get time step (should be 600 seconds from config)
dt = results['time'][1] - results['time'][0]
n_steps = len(results['time'])

print("=" * 80)
print("MASS CONSERVATION CHECK")
print("=" * 80)
print(f"Number of time steps: {n_steps}")
print(f"Time step: {dt} seconds")
print()

# Calculate initial volumes
initial_dune_x = results['dune_x'][0]
initial_dune_z = results['dune_z'][0]

# Filter out None values
initial_x_list = [initial_dune_x[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                                'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                                'ramp_toe_temp', 'ramp_toe', 'shoreline'] 
                  if k in initial_dune_x and initial_dune_x[k] is not None]
initial_z_list = [initial_dune_z[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                                'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                                'ramp_toe_temp', 'ramp_toe', 'shoreline']
                  if k in initial_dune_z and initial_dune_z[k] is not None]

initial_dune_vol = total_volume(initial_x_list, initial_z_list, substrate_x, substrate_z)
initial_vbeach = results['vbeach'][0]
initial_volbar = results['volbar'][0]

print("INITIAL STATE:")
print(f"  Dune volume: {initial_dune_vol:.4f} m³/m")
print(f"  Beach volume: {initial_vbeach:.4f} m³/m")
print(f"  Bar volume: {initial_volbar:.4f} m³/m")
print(f"  Total: {initial_dune_vol + initial_vbeach + initial_volbar:.4f} m³/m")
print()

# Track mass balance at each timestep
cumulative_error = 0.0
max_error_step = 0
max_error = 0.0

for i in range(1, n_steps):
    # Get current volumes
    curr_dune_x = results['dune_x'][i]
    curr_dune_z = results['dune_z'][i]
    
    curr_x_list = [curr_dune_x[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                             'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                             'ramp_toe_temp', 'ramp_toe', 'shoreline'] 
                   if k in curr_dune_x and curr_dune_x[k] is not None]
    curr_z_list = [curr_dune_z[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                             'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                             'ramp_toe_temp', 'ramp_toe', 'shoreline']
                   if k in curr_dune_z and curr_dune_z[k] is not None]
    
    curr_dune_vol = total_volume(curr_x_list, curr_z_list, substrate_x, substrate_z)
    
    # Get previous volumes
    prev_dune_x = results['dune_x'][i-1]
    prev_dune_z = results['dune_z'][i-1]
    
    prev_x_list = [prev_dune_x[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                             'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                             'ramp_toe_temp', 'ramp_toe', 'shoreline'] 
                   if k in prev_dune_x and prev_dune_x[k] is not None]
    prev_z_list = [prev_dune_z[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                             'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                             'ramp_toe_temp', 'ramp_toe', 'shoreline']
                   if k in prev_dune_z and prev_dune_z[k] is not None]
    
    prev_dune_vol = total_volume(prev_x_list, prev_z_list, substrate_x, substrate_z)
    
    # Calculate changes
    delta_dune = curr_dune_vol - prev_dune_vol
    delta_beach = results['vbeach'][i] - results['vbeach'][i-1]
    delta_bar = results['volbar'][i] - results['volbar'][i-1]
    
    # Calculate fluxes (these are already in m³/m, not rates)
    qs_ramp = results['qs_ramp'][i-1]  # Aeolian deposition on ramp
    qs_crest = results['qs_crest'][i-1]  # Aeolian deposition on crest
    qs_lee = results['qs_lee'][i-1]  # Aeolian deposition on lee
    qwnet = results['qwnet'][i-1]  # Net aeolian flux from beach to dune
    qbar = results['qbar'][i-1]  # Bar-beach exchange
    eroded = results['eroded_volume'][i-1]  # Eroded from dune
    
    # Expected changes
    expected_dune_change = qs_ramp + qs_crest + qs_lee - eroded
    expected_beach_change = -qwnet + qbar - (qs_ramp + qs_crest + qs_lee) + eroded
    expected_bar_change = -qbar
    
    # Mass balance error
    dune_error = delta_dune - expected_dune_change
    beach_error = delta_beach - expected_beach_change
    bar_error = delta_bar - expected_bar_change
    total_error = dune_error + beach_error + bar_error
    
    cumulative_error += total_error
    
    if abs(total_error) > abs(max_error):
        max_error = total_error
        max_error_step = i
    
    # Print details for problematic steps
    if abs(total_error) > 0.01:  # threshold
        print(f"Step {i}: Large mass balance error = {total_error:.6f} m³/m")
        print(f"  Dune change: actual={delta_dune:.6f}, expected={expected_dune_change:.6f}, error={dune_error:.6f}")
        print(f"  Beach change: actual={delta_beach:.6f}, expected={expected_beach_change:.6f}, error={beach_error:.6f}")
        print(f"  Fluxes: qs_ramp={qs_ramp:.6f}, qs_crest={qs_crest:.6f}, qs_lee={qs_lee:.6f}")
        print(f"  qwnet={qwnet:.6f}, eroded={eroded:.6f}, qbar={qbar:.6f}")
        print()

# Final volumes
final_dune_x = results['dune_x'][-1]
final_dune_z = results['dune_z'][-1]

final_x_list = [final_dune_x[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                           'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                           'ramp_toe_temp', 'ramp_toe', 'shoreline'] 
                if k in final_dune_x and final_dune_x[k] is not None]
final_z_list = [final_dune_z[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                           'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                           'ramp_toe_temp', 'ramp_toe', 'shoreline']
                if k in final_dune_z and final_dune_z[k] is not None]

final_dune_vol = total_volume(final_x_list, final_z_list, substrate_x, substrate_z)
final_vbeach = results['vbeach'][-1]
final_volbar = results['volbar'][-1]

print()
print("FINAL STATE:")
print(f"  Dune volume: {final_dune_vol:.4f} m³/m")
print(f"  Beach volume: {final_vbeach:.4f} m³/m")
print(f"  Bar volume: {final_volbar:.4f} m³/m")
print(f"  Total: {final_dune_vol + final_vbeach + final_volbar:.4f} m³/m")
print()

print("OVERALL MASS BALANCE:")
total_initial = initial_dune_vol + initial_vbeach + initial_volbar
total_final = final_dune_vol + final_vbeach + final_volbar
net_change = total_final - total_initial
print(f"  Net volume change: {net_change:.6f} m³/m")
print(f"  Cumulative mass balance error: {cumulative_error:.6f} m³/m")
print(f"  Maximum error at step {max_error_step}: {max_error:.6f} m³/m")
print(f"  Error as % of initial volume: {100 * abs(cumulative_error) / total_initial:.4f}%")
print()

if abs(cumulative_error) > 0.1:
    print("⚠️  WARNING: Significant mass balance error detected!")
else:
    print("✓ Mass is approximately conserved.")
