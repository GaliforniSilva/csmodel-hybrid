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


def calculate_compartment_volumes(dune_x, dune_z, substrate_x, substrate_z, vbeach, volbar, dclos, dfoot):
    """
    Calculate volumes for all compartments at a single timestep.
    
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
    vbeach : float
        Beach volume parameter
    volbar : float
        Bar volume
    dclos : float
        Closure depth
    dfoot : float
        Dune foot elevation
        
    Returns
    -------
    dict
        Dictionary with 'dune', 'beach', 'bar', and 'total' volumes
    """
    dune_vol = calc_dune_vol(dune_x, dune_z, substrate_x, substrate_z)
    beach_vol = calculate_actual_beach_volume(dune_x, vbeach, dclos, dfoot)
    bar_vol = volbar
    total_vol = dune_vol + beach_vol + bar_vol
    
    return {
        'dune': dune_vol,
        'beach': beach_vol,
        'bar': bar_vol,
        'total': total_vol
    }


def calculate_timestep_conservation(volumes_current, volumes_prev, fluxes, dt, initial_volumes):
    """
    Calculate conservation metrics for a single timestep.
    
    Parameters
    ----------
    volumes_current : dict
        Current timestep volumes (dune, beach, bar, total)
    volumes_prev : dict
        Previous timestep volumes (dune, beach, bar, total)
    fluxes : dict
        Fluxes for this timestep with keys:
        - 'qs_ramp': Aeolian deposition on ramp
        - 'qs_crest': Aeolian deposition on crest
        - 'qs_lee': Aeolian deposition on lee
        - 'qwnet': Net aeolian transport from beach
        - 'qbar': Bar-beach exchange
        - 'eroded_volume': Wave erosion from dune
    dt : float
        Time step in seconds
    initial_volumes : dict
        Initial volumes (for percentage calculations)
        
    Returns
    -------
    dict
        Dictionary with conservation metrics for each compartment
    """
    # Calculate actual changes
    delta_dune = volumes_current['dune'] - volumes_prev['dune']
    delta_beach = volumes_current['beach'] - volumes_prev['beach']
    delta_bar = volumes_current['bar'] - volumes_prev['bar']
    delta_total = volumes_current['total'] - volumes_prev['total']
    
    # Calculate expected changes based on fluxes
    total_aeolian = fluxes['qs_ramp'] + fluxes['qs_crest'] + fluxes['qs_lee']
    qbar_vol = fluxes['qbar'] * dt
    eroded = fluxes['eroded_volume']
    qwnet_vol = fluxes['qwnet'] * dt
    
    expected_dune_change = total_aeolian - eroded
    expected_beach_change = -qwnet_vol - qbar_vol + eroded
    expected_bar_change = qbar_vol
    expected_total_change = 0.0  # No external sources/sinks
    
    # Calculate errors
    dune_error = delta_dune - expected_dune_change
    beach_error = delta_beach - expected_beach_change
    bar_error = delta_bar - expected_bar_change
    total_error = delta_total - expected_total_change
    
    # Helper function to calculate percentage error safely
    def calc_pct_error(error, initial_value):
        return 100 * abs(error) / initial_value if initial_value > 0 else 0
    
    # Calculate percentage errors relative to initial volumes
    dune_error_pct = calc_pct_error(dune_error, initial_volumes['dune'])
    beach_error_pct = calc_pct_error(beach_error, initial_volumes['beach'])
    bar_error_pct = calc_pct_error(bar_error, initial_volumes['bar'])
    total_error_pct = calc_pct_error(total_error, initial_volumes['total'])
    
    return {
        'dune_error': dune_error,
        'beach_error': beach_error,
        'bar_error': bar_error,
        'total_error': total_error,
        'dune_error_pct': dune_error_pct,
        'beach_error_pct': beach_error_pct,
        'bar_error_pct': bar_error_pct,
        'total_error_pct': total_error_pct,
        'dune_volume': volumes_current['dune'],
        'beach_volume': volumes_current['beach'],
        'bar_volume': volumes_current['bar'],
        'total_volume': volumes_current['total']
    }


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
        - 'report_text': full formatted report as a string
    """
    
    # Build report text
    report_lines = []
    
    def report_print(text=""):
        """Helper function to both print and append to report"""
        report_lines.append(text)
        if verbose:
            print(text)
    
    report_print("=" * 80)
    report_print("CS MODEL MASS CONSERVATION REPORT")
    report_print("=" * 80)
    report_print()
    
    # Get parameters from config
    dclos = config['erosion']['dclos']
    dfoot = config['dune_parameters']['dfoot']
    
    substrate_x = results['substrate_x']
    substrate_z = results['substrate_z']
    n_steps = len(results['time'])
    dt = results['time'][1] - results['time'][0] if n_steps > 1 else 0
    total_time_hours = results['time'][-1] / 3600 if n_steps > 0 else 0
    
    report_print(f"Simulation details:")
    report_print(f"  Number of time steps: {n_steps}")
    report_print(f"  Time step: {dt:.0f} seconds")
    report_print(f"  Total simulation time: {total_time_hours:.1f} hours ({total_time_hours/24:.1f} days)")
    report_print()
    
    # Calculate volumes at start and end
    initial_dune_vol = calc_dune_vol(results['dune_x'][0], results['dune_z'][0], substrate_x, substrate_z)
    initial_beach_vol = calculate_actual_beach_volume(results['dune_x'][0], results['vbeach'][0], dclos, dfoot)
    initial_bar_vol = results['volbar'][0]
    initial_total = initial_dune_vol + initial_beach_vol + initial_bar_vol
    
    final_dune_vol = calc_dune_vol(results['dune_x'][-1], results['dune_z'][-1], substrate_x, substrate_z)
    final_beach_vol = calculate_actual_beach_volume(results['dune_x'][-1], results['vbeach'][-1], dclos, dfoot)
    final_bar_vol = results['volbar'][-1]
    final_total = final_dune_vol + final_beach_vol + final_bar_vol
    
    report_print("INITIAL STATE:")
    report_print(f"  Dune volume:  {initial_dune_vol:>12.4f} m³/m")
    report_print(f"  Beach volume: {initial_beach_vol:>12.4f} m³/m")
    report_print(f"  Bar volume:   {initial_bar_vol:>12.4f} m³/m")
    report_print(f"  Total volume: {initial_total:>12.4f} m³/m")
    report_print()

    report_print("FINAL STATE:")
    report_print(f"  Dune volume:  {final_dune_vol:>12.4f} m³/m")
    report_print(f"  Beach volume: {final_beach_vol:>12.4f} m³/m")
    report_print(f"  Bar volume:   {final_bar_vol:>12.4f} m³/m")
    report_print(f"  Total volume: {final_total:>12.4f} m³/m")
    report_print()
    
    # Calculate changes
    delta_dune = final_dune_vol - initial_dune_vol
    delta_beach = final_beach_vol - initial_beach_vol
    delta_bar = final_bar_vol - initial_bar_vol
    delta_total = final_total - initial_total
    
    report_print("VOLUME CHANGES:")
    report_print(f"  Δ Dune:  {delta_dune:>+12.4f} m³/m ({100*delta_dune/initial_dune_vol:>+7.2f}%)")
    report_print(f"  Δ Beach: {delta_beach:>+12.4f} m³/m ({100*delta_beach/initial_beach_vol:>+7.2f}%)")
    report_print(f"  Δ Bar:   {delta_bar:>+12.4f} m³/m ({100*delta_bar/initial_bar_vol:>+7.2f}%)")
    report_print(f"  Δ Total: {delta_total:>+12.4f} m³/m ({100*delta_total/initial_total:>+7.2f}%)")
    report_print()
    
    # Calculate integrated fluxes
    total_qs_ramp = np.sum(results['qs_ramp'][:n_steps])
    total_qs_crest = np.sum(results['qs_crest'][:n_steps])
    total_qs_lee = np.sum(results['qs_lee'][:n_steps])
    total_aeolian = total_qs_ramp + total_qs_crest + total_qs_lee
    
    total_qwnet_vol = np.sum([results['qwnet'][i] * dt for i in range(n_steps)])
    total_qbar_vol = np.sum([results['qbar'][i] * dt for i in range(n_steps)])
    total_eroded = np.sum(results['eroded_volume'][:n_steps])
    
    report_print("INTEGRATED FLUXES (over entire simulation):")
    report_print(f"  Aeolian transport from beach to dune: {total_aeolian:>12.4f} m³/m")
    report_print(f"    - Deposited on ramp:  {total_qs_ramp:>12.4f} m³/m")
    report_print(f"    - Deposited on crest: {total_qs_crest:>12.4f} m³/m")
    report_print(f"    - Deposited on lee:   {total_qs_lee:>12.4f} m³/m")
    report_print(f"  Bar-beach exchange (to bar if +):     {total_qbar_vol:>12.4f} m³/m")
    report_print(f"  Wave erosion (from dune):             {-total_eroded:>12.4f} m³/m")
    report_print()
    
    # Mass balance assessment
    expected_dune_change = total_aeolian - total_eroded
    expected_beach_change = -total_aeolian - total_qbar_vol + total_eroded
    expected_bar_change = total_qbar_vol
    expected_total_change = 0.0  # No external sources/sinks
    
    report_print("MASS BALANCE ASSESSMENT:")
    report_print()
    report_print("Expected changes (based on fluxes):")
    report_print(f"  Expected Δ Dune:  {expected_dune_change:>+12.4f} m³/m")
    report_print(f"  Expected Δ Beach: {expected_beach_change:>+12.4f} m³/m")
    report_print(f"  Expected Δ Bar:   {expected_bar_change:>+12.4f} m³/m")
    report_print(f"  Expected Δ Total: {expected_total_change:>+12.4f} m³/m")
    report_print()
    
    dune_error = delta_dune - expected_dune_change
    beach_error = delta_beach - expected_beach_change
    bar_error = delta_bar - expected_bar_change
    total_error = delta_total - expected_total_change
    
    report_print("Mass balance errors:")
    report_print(f"  Dune error:  {dune_error:>+12.6f} m³/m ({100*abs(dune_error)/initial_dune_vol:.4f}%)")
    report_print(f"  Beach error: {beach_error:>+12.6f} m³/m ({100*abs(beach_error)/initial_beach_vol:.4f}%)")
    report_print(f"  Bar error:   {bar_error:>+12.6f} m³/m ({100*abs(bar_error)/initial_bar_vol:.4f}%)")
    report_print(f"  TOTAL ERROR: {total_error:>+12.6f} m³/m ({100*abs(total_error)/initial_total:.6f}%)")
    report_print()
    
    # Time-series conservation analysis (if available)
    if 'conservation' in results and len(results['conservation'].get('total_error', [])) > 0:
        report_print("TIME-SERIES CONSERVATION ANALYSIS:")
        report_print()
        
        conservation = results['conservation']
        
        # Find maximum errors
        max_dune_error = np.max(np.abs(conservation['dune_error']))
        max_beach_error = np.max(np.abs(conservation['beach_error']))
        max_bar_error = np.max(np.abs(conservation['bar_error']))
        max_total_error = np.max(np.abs(conservation['total_error']))
        
        # Find when maximum errors occurred
        idx_max_dune = np.argmax(np.abs(conservation['dune_error']))
        idx_max_beach = np.argmax(np.abs(conservation['beach_error']))
        idx_max_bar = np.argmax(np.abs(conservation['bar_error']))
        idx_max_total = np.argmax(np.abs(conservation['total_error']))
        
        # Calculate time for each max error
        time_max_dune = idx_max_dune * dt / 3600
        time_max_beach = idx_max_beach * dt / 3600
        time_max_bar = idx_max_bar * dt / 3600
        time_max_total = idx_max_total * dt / 3600
        
        report_print("Maximum errors during simulation:")
        report_print(f"  Dune:  {max_dune_error:>12.6e} m³/m at t={time_max_dune:>7.2f}h (step {idx_max_dune})")
        report_print(f"  Beach: {max_beach_error:>12.6e} m³/m at t={time_max_beach:>7.2f}h (step {idx_max_beach})")
        report_print(f"  Bar:   {max_bar_error:>12.6e} m³/m at t={time_max_bar:>7.2f}h (step {idx_max_bar})")
        report_print(f"  Total: {max_total_error:>12.6e} m³/m at t={time_max_total:>7.2f}h (step {idx_max_total})")
        report_print()
        
        # Calculate cumulative absolute errors
        cumulative_dune = np.sum(np.abs(conservation['dune_error']))
        cumulative_beach = np.sum(np.abs(conservation['beach_error']))
        cumulative_bar = np.sum(np.abs(conservation['bar_error']))
        cumulative_total = np.sum(np.abs(conservation['total_error']))
        
        report_print("Cumulative absolute errors:")
        report_print(f"  Dune:  {cumulative_dune:>12.6e} m³/m")
        report_print(f"  Beach: {cumulative_beach:>12.6e} m³/m")
        report_print(f"  Bar:   {cumulative_bar:>12.6e} m³/m")
        report_print(f"  Total: {cumulative_total:>12.6e} m³/m")
        report_print()
        
        # Identify problematic timesteps (errors > 0.01% of initial volume)
        threshold_pct = 0.01
        problematic_steps = []
        for i in range(len(conservation['total_error'])):
            if (conservation['dune_error_pct'][i] > threshold_pct or 
                conservation['beach_error_pct'][i] > threshold_pct or 
                conservation['bar_error_pct'][i] > threshold_pct):
                problematic_steps.append(i)
        
        if problematic_steps:
            n_show = min(5, len(problematic_steps))
            report_print(f"Timesteps with errors > {threshold_pct}% (showing first {n_show} of {len(problematic_steps)}):")
            for step in problematic_steps[:n_show]:
                time_h = step * dt / 3600
                report_print(f"  Step {step:4d} (t={time_h:7.2f}h): "
                          f"Dune={conservation['dune_error_pct'][step]:.4f}%, "
                          f"Beach={conservation['beach_error_pct'][step]:.4f}%, "
                          f"Bar={conservation['bar_error_pct'][step]:.4f}%")
        else:
            report_print(f"✓ No timesteps with errors > {threshold_pct}% of initial volumes")
        
        report_print()
    
    # Assessment
    error_threshold = 0.1  # 0.1% threshold
    total_error_percent = abs(100*total_error/initial_total)
    conserved = total_error_percent < error_threshold
    
    if conserved:
        report_print("✓ RESULT: Mass is conserved within acceptable limits!")
        report_print(f"  (Error < {error_threshold}% threshold)")
    else:
        report_print("⚠ WARNING: Mass balance error exceeds acceptable threshold!")
        report_print(f"  (Error > {error_threshold}% threshold)")
    
    report_print()
    report_print("=" * 80)
    
    # Join all report lines into a single string
    report_text = "\n".join(report_lines)
    
    return {
        'conserved': conserved,
        'total_error_percent': total_error_percent,
        'total_error': total_error,
        'initial_total': initial_total,
        'final_total': final_total,
        'dune_error': dune_error,
        'beach_error': beach_error,
        'bar_error': bar_error,
        'error_threshold': error_threshold,
        'report_text': report_text
    }
