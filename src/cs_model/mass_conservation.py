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

from cs_model.functions import partial_volume


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

    The volume is integrated from ``lee_foot`` to ``ramp_toe`` using
    ``partial_volume``.  Starting the integration at ``lee_foot`` (rather than
    at ``back_ref``) avoids a phantom-volume artefact: when the lee foot rides
    up the dike back face (``lee_foot.z > dfoot``), the straight line drawn by
    ``total_volume`` from ``back_ref`` (at z = dfoot) to the elevated
    ``lee_foot`` passes above the flat back-shore substrate (also at z = dfoot),
    creating large spurious volume that is not backed by any real sediment flux.
    Because the back-shore from ``back_ref`` to the base of the lee slope is
    always at z = dfoot = substrate elevation, it contributes zero to the true
    dune volume; ``partial_volume(lee_foot, ramp_toe)`` computes the same
    correct value in all cases.

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
        Total dune volume above substrate (m³/m)
    """
    return partial_volume(dune_x, dune_z, substrate_x, substrate_z,
                          dune_x['lee_foot'], dune_x['ramp_toe'])


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


def check_mass_conservation(results, config, verbose=True, skip_first_step=False):
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
    skip_first_step : bool, optional
        If True, the first model time step (i=0) is excluded from the
        conservation analysis.  The initial state is then taken as the
        state *after* the first step (``results['dune_x'][1]``,
        ``results['vbeach'][0]``, ``results['volbar'][0]``) and all
        flux integrals are accumulated starting from step index 1.
        Default is False (use the true initial state, i.e. current
        behaviour).
        
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

    report_lines = []

    def report_print(text=""):
        report_lines.append(text)
        if verbose:
            print(text)

    report_print("=" * 80)
    report_print("CS MODEL MASS CONSERVATION REPORT")
    report_print("=" * 80)
    report_print()

    # ------------------------------------------------------------------
    # PARAMETERS
    # ------------------------------------------------------------------
    dclos = config['erosion']['dclos']
    dfoot = config['dune_parameters']['dfoot']
    al    = config['constants']['al']
    bl    = config['constants']['bl']

    substrate_x = results['substrate_x']
    substrate_z = results['substrate_z']

    n_steps = len(results['time'])
    dt = results['time'][1] - results['time'][0] if n_steps > 1 else 0.0
    total_time_hours = results['time'][-1] / 3600 if n_steps > 0 else 0.0

    report_print(f"Simulation details:")
    report_print(f"  Number of time steps: {n_steps}")
    report_print(f"  Time step: {dt:.0f} s")
    report_print(f"  Total simulation time: {total_time_hours:.1f} h "
                 f"({total_time_hours/24:.1f} days)")
    if skip_first_step:
        report_print("  (First timestep excluded from conservation analysis)")
    report_print()

    flux_start = 1 if skip_first_step else 0

    # ------------------------------------------------------------------
    # INITIAL STATE
    # ------------------------------------------------------------------
    if skip_first_step:
        initial_vbeach = results['vbeach'][0]
        initial_dune_vol = calc_dune_vol(
            results['dune_x'][1], results['dune_z'][1],
            substrate_x, substrate_z
        )
        initial_beach_vol = calculate_actual_beach_volume(
            results['dune_x'][1], initial_vbeach, dclos, dfoot
        )
        initial_bar_vol = results['volbar'][0]
    else:
        initial_vbeach = (
            (results['dune_x'][0]['shoreline']
             - results['dune_x'][0]['ramp_toe']) - bl
        ) / al
        initial_dune_vol = calc_dune_vol(
            results['dune_x'][0], results['dune_z'][0],
            substrate_x, substrate_z
        )
        initial_beach_vol = calculate_actual_beach_volume(
            results['dune_x'][0], initial_vbeach, dclos, dfoot
        )
        initial_bar_vol = (
            results['volbar'][0] - results['qbar'][0] * dt
        )

    initial_total = (
        initial_dune_vol + initial_beach_vol + initial_bar_vol
    )

    # ------------------------------------------------------------------
    # FINAL STATE
    # ------------------------------------------------------------------
    final_dune_vol = calc_dune_vol(
        results['dune_x'][-1], results['dune_z'][-1],
        substrate_x, substrate_z
    )
    final_beach_vol = calculate_actual_beach_volume(
        results['dune_x'][-1], results['vbeach'][-1], dclos, dfoot
    )
    final_bar_vol = results['volbar'][-1]
    final_total = final_dune_vol + final_beach_vol + final_bar_vol

    report_print("INITIAL TOTAL VOLUME:")
    report_print(f"  {initial_total:.6f} m³/m")
    report_print("FINAL TOTAL VOLUME:")
    report_print(f"  {final_total:.6f} m³/m")
    report_print()

    delta_total = final_total - initial_total

    # ------------------------------------------------------------------
    # INTEGRATED FLUXES
    # ------------------------------------------------------------------
    total_qs_ramp  = np.sum(results['qs_ramp'][flux_start:n_steps])
    total_qs_crest = np.sum(results['qs_crest'][flux_start:n_steps])
    total_qs_lee   = np.sum(results['qs_lee'][flux_start:n_steps])
    total_aeolian  = total_qs_ramp + total_qs_crest + total_qs_lee

    total_qbar_vol = np.sum(np.asarray(results['qbar'][flux_start:n_steps]) * dt)
    total_eroded   = np.sum(results['eroded_volume'][flux_start:n_steps])

    total_transport = (
        abs(total_aeolian) +
        abs(total_qbar_vol) +
        abs(total_eroded)
    )

    # ------------------------------------------------------------------
    # PHANTOM BEACH CORRECTION
    # ------------------------------------------------------------------
    phantom_correction = 0.0
    vb_prev = initial_vbeach

    for i in range(flux_start, n_steps):
        rt_old = results['dune_x'][i]['ramp_toe']
        sl_old = results['dune_x'][i]['shoreline']
        rt_new = results['dune_x'][i + 1]['ramp_toe']

        beach_width_old = sl_old - rt_old
        db_old = vb_prev / beach_width_old if beach_width_old > 0 else 0.0

        phantom_correction += (rt_new - rt_old) * (dfoot - db_old)
        vb_prev = results['vbeach'][i]

    # ------------------------------------------------------------------
    # EXPECTED TOTAL CHANGE
    # ------------------------------------------------------------------
    expected_total_change = phantom_correction

    total_error = delta_total - expected_total_change
    total_error_percent = abs(100 * total_error / initial_total)
    error_transport_rel = (
        abs(total_error) / total_transport
        if total_transport > 0 else 0.0
    )

    drift_per_step = total_error / n_steps
    drift_per_hour = total_error / total_time_hours if total_time_hours > 0 else 0.0

    # Per-compartment breakdown
    delta_dune  = final_dune_vol  - initial_dune_vol
    delta_beach = final_beach_vol - initial_beach_vol
    delta_bar   = final_bar_vol   - initial_bar_vol

    exp_dune  = total_aeolian - total_eroded
    # Beach: loses aeolian flux to dune (-total_aeolian), receives eroded dune material
    # (+total_eroded), loses sediment to bar (-total_qbar_vol), and changes by a
    # geometric phantom when the dune toe (ramp_toe) shifts (+phantom_correction).
    exp_beach = -total_aeolian - total_qbar_vol + total_eroded + phantom_correction
    exp_bar   = total_qbar_vol

    dune_error  = delta_dune  - exp_dune
    beach_error = delta_beach - exp_beach
    bar_error   = delta_bar   - exp_bar

    def _pct(err, ref):
        """Return |err| as a percentage of ref, or 0 if ref is non-positive."""
        return 100 * abs(err) / ref if ref > 0 else 0.0

    # ------------------------------------------------------------------
    # REPORT
    # ------------------------------------------------------------------
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
    report_print("VOLUME CHANGES:")
    report_print(f"  Δ Dune:  {delta_dune:>+12.4f} m³/m  "
                 f"(expected {exp_dune:>+.4f}  error {dune_error:>+.4f}  "
                 f"{_pct(dune_error, initial_dune_vol):.3f}%)")
    report_print(f"  Δ Beach: {delta_beach:>+12.4f} m³/m  "
                 f"(expected {exp_beach:>+.4f}  error {beach_error:>+.4f}  "
                 f"{_pct(beach_error, initial_beach_vol):.3f}%)")
    report_print(f"  Δ Bar:   {delta_bar:>+12.4f} m³/m  "
                 f"(expected {exp_bar:>+.4f}  error {bar_error:>+.4f}  "
                 f"{_pct(bar_error, initial_bar_vol):.3f}%)")
    report_print(f"  Δ Total: {delta_total:>+12.4f} m³/m  "
                 f"(phantom correction {phantom_correction:>+.4f})")
    report_print()
    report_print("INTEGRATED FLUXES:")
    report_print(f"  Aeolian (ramp+crest+lee): {total_aeolian:>+12.4f} m³/m")
    report_print(f"    - ramp:  {total_qs_ramp:>+12.4f} m³/m")
    report_print(f"    - crest: {total_qs_crest:>+12.4f} m³/m")
    report_print(f"    - lee:   {total_qs_lee:>+12.4f} m³/m")
    report_print(f"  Bar-beach exchange (+ = to bar): {total_qbar_vol:>+12.4f} m³/m")
    report_print(f"  Wave erosion from dune:          {total_eroded:>+12.4f} m³/m")
    report_print()
    report_print("MASS BALANCE ASSESSMENT:")
    report_print()
    report_print(f"  Total mass error: {total_error:+.6e} m³/m")
    report_print(f"  Relative to initial volume: {total_error_percent:.4f}%")
    report_print(f"  Relative to total transport: {100*error_transport_rel:.4f}%")
    report_print()
    report_print("DRIFT DIAGNOSTICS (informational):")
    report_print(f"  Mean drift per step: {drift_per_step:+.3e} m³/m")
    report_print(f"  Mean drift rate:     {drift_per_hour:+.3e} m³/m per hour")
    report_print()

    # ------------------------------------------------------------------
    # VERDICT
    # ------------------------------------------------------------------
    error_threshold = 0.1  # percent of initial volume
    conserved = total_error_percent < error_threshold

    if conserved:
        report_print("✓ RESULT: Mass is conserved within acceptable limits.")
        report_print(f"  (|ΔV| < {error_threshold}% of initial total volume)")
    else:
        report_print("⚠ WARNING: Mass balance error exceeds acceptable limits.")
        report_print(f"  (|ΔV| > {error_threshold}% threshold)")

    report_print()
    report_print("=" * 80)

    return {
        'conserved': conserved,
        'total_error_percent': total_error_percent,
        'total_error': total_error,
        'dune_error': dune_error,
        'beach_error': beach_error,
        'bar_error': bar_error,
        'transport_relative_error': error_transport_rel,
        'drift_per_step': drift_per_step,
        'drift_per_hour': drift_per_hour,
        'initial_total': initial_total,
        'final_total': final_total,
        'error_threshold': error_threshold,
        'report_text': "\n".join(report_lines)
    }

