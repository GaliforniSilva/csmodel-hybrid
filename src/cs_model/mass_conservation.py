"""
mass_conservation.py
Authors:
    Filipe Galiforni Silva (TU Delft)

TU Delft, 2026

Description:
    Mass conservation analysis for CS Model simulations.
    This module provides functions to check and report on mass conservation
    in simulation results.

License:
    Copyright (c) 2026, Delft University of Technology.
    All rights reserved.
"""

import numpy as np
import sys
import os

# Import from local package
package_dir = os.path.dirname(os.path.abspath(__file__))
if package_dir not in sys.path:
    sys.path.insert(0, package_dir)

from cs_model.functions import total_volume


def calculate_actual_beach_volume(dune_x, vbeach, dclos, dfoot):
    """
    Calculate the actual beach volume (vbeachtot) from geometry.
    
    Parameters
    ----------
    dune_x : dict
        Dictionary of dune x-coordinates
    vbeach : float
        Beach volume parameter
    dclos : float
        Closure depth
    dfoot : float
        Dune foot elevation
        
    Returns
    -------
    float
        Actual beach volume
    """
    shoreline = dune_x['shoreline']
    ramp_toe = dune_x['ramp_toe']
    db = vbeach / (shoreline - ramp_toe)
    vbeachtot = ramp_toe * (dclos + dfoot) + (shoreline - ramp_toe) * (dclos + db)
    return vbeachtot


def calc_dune_vol(dune_x, dune_z, substrate_x, substrate_z):
    """
    Calculate total dune volume above substrate.
    
    Parameters
    ----------
    dune_x : dict
        Dictionary of dune x-coordinates
    dune_z : dict
        Dictionary of dune z-coordinates
    substrate_x : list
        List of substrate x-coordinates
    substrate_z : list
        List of substrate z-coordinates
        
    Returns
    -------
    float
        Total dune volume
    """
    x_list = [dune_x[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                   'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                   'ramp_toe_temp', 'ramp_toe', 'shoreline'] 
              if k in dune_x and dune_x[k] is not None]
    z_list = [dune_z[k] for k in ['back_ref', 'lee_foot', 'lee_crest', 'crest_temp', 
                                   'seaward_crest', 'ramp_dike_temp', 'ramp_crest', 
                                   'ramp_toe_temp', 'ramp_toe', 'shoreline']
              if k in dune_z and dune_z[k] is not None]
    return total_volume(x_list, z_list, substrate_x, substrate_z)


def check_mass_conservation(results, config, verbose=True):
    """
    Check mass conservation in model results and generate a report.
    
    Parameters
    ----------
    results : dict
        Model output dictionary containing time series data
    config : dict
        Configuration dictionary with model parameters
    verbose : bool, optional
        If True, print detailed report. Default is True.
        
    Returns
    -------
    dict
        Dictionary containing mass conservation metrics including:
        - 'conserved': bool indicating if mass is conserved within threshold
        - 'total_error_percent': total mass balance error as percentage
        - 'total_error': absolute total mass balance error
        - 'initial_total': initial total volume
        - 'final_total': final total volume
    """
    
    if verbose:
        print("=" * 80)
        print("CS MODEL MASS CONSERVATION REPORT")
        print("=" * 80)
        print()
    
    # Get parameters from config
    dclos = config['erosion']['dclos']
    dfoot = config['dune_parameters']['dfoot']
    
    substrate_x = results['substrate_x']
    substrate_z = results['substrate_z']
    n_steps = len(results['time'])
    dt = results['time'][1] - results['time'][0] if n_steps > 1 else 0
    total_time_hours = results['time'][-1] / 3600 if n_steps > 0 else 0
    
    if verbose:
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
    
    if verbose:
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
    
    if verbose:
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
    
    if verbose:
        print("INTEGRATED FLUXES (over entire simulation):")
        print(f"  Aeolian transport from beach to dune: {total_aeolian:>12.4f} m³/m")
        print(f"    - Deposited on ramp:  {total_qs_ramp:>12.4f} m³/m")
        print(f"    - Deposited on crest: {total_qs_crest:>12.4f} m³/m")
        print(f"    - Deposited on lee:   {total_qs_lee:>12.4f} m³/m")
        print(f"  Bar-beach exchange (to bar if +):     {total_qbar_vol:>12.4f} m³/m")
        print(f"  Wave erosion (from dune):             {-total_eroded:>12.4f} m³/m")
        print()
    
    # Mass balance assessment
    expected_dune_change = total_aeolian - total_eroded
    expected_beach_change = -total_aeolian - total_qbar_vol + total_eroded
    expected_bar_change = total_qbar_vol
    expected_total_change = 0.0  # No external sources/sinks
    
    if verbose:
        print("MASS BALANCE ASSESSMENT:")
        print()
        print("Expected changes (based on fluxes):")
        print(f"  Expected Δ Dune:  {expected_dune_change:>+12.4f} m³/m")
        print(f"  Expected Δ Beach: {expected_beach_change:>+12.4f} m³/m")
        print(f"  Expected Δ Bar:   {expected_bar_change:>+12.4f} m³/m")
        print(f"  Expected Δ Total: {expected_total_change:>+12.4f} m³/m")
        print()
    
    dune_error = delta_dune - expected_dune_change
    beach_error = delta_beach - expected_beach_change
    bar_error = delta_bar - expected_bar_change
    total_error = delta_total - expected_total_change
    
    if verbose:
        print("Mass balance errors:")
        print(f"  Dune error:  {dune_error:>+12.6f} m³/m ({100*abs(dune_error)/initial_dune_vol:.4f}%)")
        print(f"  Beach error: {beach_error:>+12.6f} m³/m ({100*abs(beach_error)/initial_beach_vol:.4f}%)")
        print(f"  Bar error:   {bar_error:>+12.6f} m³/m ({100*abs(bar_error)/initial_bar_vol:.4f}%)")
        print(f"  TOTAL ERROR: {total_error:>+12.6f} m³/m ({100*abs(total_error)/initial_total:.6f}%)")
        print()
    
    # Assessment
    error_threshold = 0.1  # 0.1% threshold
    total_error_percent = abs(100*total_error/initial_total)
    conserved = total_error_percent < error_threshold
    
    if verbose:
        if conserved:
            print("✓ RESULT: Mass is conserved within acceptable limits!")
            print(f"  (Error < {error_threshold}% threshold)")
        else:
            print("⚠ WARNING: Mass balance error exceeds acceptable threshold!")
            print(f"  (Error > {error_threshold}% threshold)")
        
        print()
        print("=" * 80)
    
    return {
        'conserved': conserved,
        'total_error_percent': total_error_percent,
        'total_error': total_error,
        'initial_total': initial_total,
        'final_total': final_total,
        'dune_error': dune_error,
        'beach_error': beach_error,
        'bar_error': bar_error,
        'error_threshold': error_threshold
    }
