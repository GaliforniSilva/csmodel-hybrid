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
        - 'nour_dune': Nourishment added to dune (optional, default 0)
        - 'nour_beach': Nourishment added to beach (optional, default 0)
        - 'nour_bar': Nourishment added to bar (optional, default 0)
    dt : float
        Time step in seconds
    initial_volumes : dict
        Initial volumes (used as fallback denominator when per-step transport is zero)
        
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
    
    # Nourishment (external sediment sources, zero if not provided)
    nour_dune  = fluxes.get('nour_dune', 0.0)
    nour_beach = fluxes.get('nour_beach', 0.0)
    nour_bar   = fluxes.get('nour_bar', 0.0)

    # Calculate expected changes based on fluxes
    total_aeolian = fluxes['qs_ramp'] + fluxes['qs_crest'] + fluxes['qs_lee']
    qbar_vol = fluxes['qbar'] * dt
    eroded = fluxes['eroded_volume']
    qwnet_vol = fluxes['qwnet'] * dt
    
    expected_dune_change  = total_aeolian - eroded + nour_dune
    expected_beach_change = -qwnet_vol - qbar_vol + eroded + nour_beach
    expected_bar_change   = qbar_vol + nour_bar
    expected_total_change = nour_dune + nour_beach + nour_bar  # internal fluxes cancel
    
    # Calculate errors
    dune_error  = delta_dune  - expected_dune_change
    beach_error = delta_beach - expected_beach_change
    bar_error   = delta_bar   - expected_bar_change
    total_error = delta_total - expected_total_change
    
    # Per-step transport magnitude: used to normalize the total error percentage so
    # the metric scales with simulation activity rather than with initial volume.
    step_transport = (
        abs(total_aeolian) + abs(qwnet_vol) + abs(qbar_vol) + abs(eroded)
        + abs(nour_dune) + abs(nour_beach) + abs(nour_bar)
    )

    # Helper function to calculate percentage error safely
    def calc_pct_error(error, denominator):
        return 100 * abs(error) / denominator if denominator > 0 else 0

    # Compartment errors: relative to initial compartment volume (diagnostic)
    dune_error_pct  = calc_pct_error(dune_error,  initial_volumes['dune'])
    beach_error_pct = calc_pct_error(beach_error, initial_volumes['beach'])
    bar_error_pct   = calc_pct_error(bar_error,   initial_volumes['bar'])

    # Total error: relative to per-step transport (scales with duration & activity)
    total_error_pct = calc_pct_error(total_error, step_transport)
    
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

    # Nourishment: external sediment added to the system — counts toward expected change
    if 'nourishment' in results:
        nour = results['nourishment']
        total_nour_dune  = np.sum(nour['dune'][flux_start:n_steps])
        total_nour_beach = np.sum(nour['beach'][flux_start:n_steps])
        total_nour_bar   = np.sum(nour['bar'][flux_start:n_steps])
        total_nourishment = total_nour_dune + total_nour_beach + total_nour_bar
    else:
        total_nourishment = 0.0

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
    # Phantom beach correction (geometric effect of ramp_toe movement) plus
    # any externally added sediment (nourishments) together define the expected
    # change in total volume.  Internal fluxes cancel exactly by design.
    expected_total_change = phantom_correction + total_nourishment

    total_error = delta_total - expected_total_change
    total_error_percent = abs(100 * total_error / initial_total)
    error_transport_rel = (
        abs(total_error) / total_transport
        if total_transport > 0 else 0.0
    )

    drift_per_step = total_error / n_steps
    drift_per_hour = total_error / total_time_hours if total_time_hours > 0 else 0.0
    drift_per_year = drift_per_hour * 8760  # annualised drift rate

    # ------------------------------------------------------------------
    # REPORT
    # ------------------------------------------------------------------
    if total_nourishment != 0.0:
        report_print(f"  Total nourishment applied: {total_nourishment:+.4f} m³/m")
        report_print()

    report_print("MASS BALANCE ASSESSMENT:")
    report_print()
    report_print(f"  Total mass error: {total_error:+.6e} m³/m")
    report_print(f"  Relative to initial volume:   {total_error_percent:.4f}%")
    report_print(f"  Relative to total transport:  {100*error_transport_rel:.4f}%")
    report_print()
    report_print("DRIFT DIAGNOSTICS (informational):")
    report_print(f"  Mean drift per step: {drift_per_step:+.3e} m³/m")
    report_print(f"  Mean drift rate:     {drift_per_hour:+.3e} m³/m per hour")
    report_print(f"  Annualised drift:    {drift_per_year:+.3e} m³/m per year")
    report_print()

    # ------------------------------------------------------------------
    # VERDICT  (dual criterion — scales with both initial volume and activity)
    # ------------------------------------------------------------------
    # Primary: error as fraction of total transport — scales correctly with
    # simulation duration and activity level.
    # Secondary: error as fraction of initial volume — absolute sanity check.
    # Both thresholds must be exceeded for a FAILED result, so that a model
    # with tiny transport and a small absolute error is not penalised, and a
    # long run with huge transport is still checked against an absolute bound.
    transport_threshold = 0.5   # % of total transport (primary)
    volume_threshold    = 1.0   # % of initial total volume (secondary)

    transport_ok = (100 * error_transport_rel) < transport_threshold
    volume_ok    = total_error_percent < volume_threshold
    conserved    = transport_ok or volume_ok  # pass if *either* criterion is met

    if conserved:
        report_print("✓ RESULT: Mass is conserved within acceptable limits.")
        if transport_ok:
            report_print(f"  (|ΔV|/transport < {transport_threshold}% — primary criterion met)")
        if volume_ok:
            report_print(f"  (|ΔV|/V₀ < {volume_threshold}% — secondary criterion met)")
    else:
        report_print("⚠ WARNING: Mass balance error exceeds acceptable limits.")
        report_print(f"  (|ΔV|/transport = {100*error_transport_rel:.4f}% ≥ {transport_threshold}% primary threshold)")
        report_print(f"  (|ΔV|/V₀ = {total_error_percent:.4f}% ≥ {volume_threshold}% secondary threshold)")

    report_print()
    report_print("=" * 80)

    # ------------------------------------------------------------------
    # TIMESTEP LOCALISATION DIAGNOSTICS (only if conservation fails)
    # ------------------------------------------------------------------
    if not conserved and 'conservation' in results:
        conservation = results['conservation']
    
        if 'total_error' in conservation and len(conservation['total_error']) > 1:
            total_err = np.asarray(conservation['total_error'])
    
            # Incremental error per timestep
            inc_err = np.diff(total_err)
    
            # Characteristic scale for "normal" error
            err_scale = np.median(np.abs(inc_err)) + 1e-14
    
            # Flag spikes: > N times typical increment
            spike_factor = 20.0
            spike_idx = np.where(np.abs(inc_err) > spike_factor * err_scale)[0]
    
            report_print("TIMESTEP DIAGNOSTICS (non-systemic errors):")
            report_print()
    
            if len(spike_idx) == 0:
                report_print("  No isolated error spikes detected.")
                report_print("  → Error likely due to systematic drift.")
            else:
                report_print(
                    f"  Detected {len(spike_idx)} timestep(s) with abnormal mass imbalance:"
                )
    
                n_show = min(10, len(spike_idx))
                for k in spike_idx[:n_show]:
                    t_h = (k + 1) * dt / 3600
                    report_print(
                        f"    Step {k+1:6d} (t = {t_h:8.2f} h): "
                        f"Δerror = {inc_err[k]:+.3e} m³/m"
                    )
    
                if len(spike_idx) > n_show:
                    report_print(
                        f"    ... ({len(spike_idx) - n_show} more not shown)"
                    )
    
            report_print()

    return {
        'conserved': conserved,
        'total_error_percent': total_error_percent,
        'total_error': total_error,
        'transport_relative_error': error_transport_rel,
        'drift_per_step': drift_per_step,
        'drift_per_hour': drift_per_hour,
        'drift_per_year': drift_per_year,
        'initial_total': initial_total,
        'final_total': final_total,
        'transport_threshold': transport_threshold,
        'volume_threshold': volume_threshold,
        'report_text': "\n".join(report_lines)
    }    
    # # Build report text
    # report_lines = []
    
    # def report_print(text=""):
    #     """Helper function to both print and append to report"""
    #     report_lines.append(text)
    #     if verbose:
    #         print(text)
    
    # report_print("=" * 80)
    # report_print("CS MODEL MASS CONSERVATION REPORT")
    # report_print("=" * 80)
    # report_print()
    
    # # Get parameters from config
    # dclos = config['erosion']['dclos']
    # dfoot = config['dune_parameters']['dfoot']
    # al    = config['constants']['al']
    # bl    = config['constants']['bl']
    
    # substrate_x = results['substrate_x']
    # substrate_z = results['substrate_z']
    # n_steps = len(results['time'])
    # dt = results['time'][1] - results['time'][0] if n_steps > 1 else 0
    # total_time_hours = results['time'][-1] / 3600 if n_steps > 0 else 0
    
    # report_print(f"Simulation details:")
    # report_print(f"  Number of time steps: {n_steps}")
    # report_print(f"  Time step: {dt:.0f} seconds")
    # report_print(f"  Total simulation time: {total_time_hours:.1f} hours ({total_time_hours/24:.1f} days)")
    # if skip_first_step:
    #     report_print(f"  (First time step excluded from conservation analysis)")
    # report_print()
    
    # # Determine the flux/phantom-correction start index.
    # # When skip_first_step=True, step 0 is excluded and the "initial state" is
    # # taken as the model state *after* the first time step.
    # flux_start = 1 if skip_first_step else 0
    
    # # Calculate volumes at start and end.
    # if skip_first_step:
    #     # Use state after first step directly (no reconstruction needed).
    #     # results['dune_x'][1] is the profile after step 0; results['vbeach'][0]
    #     # and results['volbar'][0] are the corresponding beach/bar volumes.
    #     initial_vbeach = results['vbeach'][0]
    #     initial_dune_vol = calc_dune_vol(results['dune_x'][1], results['dune_z'][1], substrate_x, substrate_z)
    #     initial_beach_vol = calculate_actual_beach_volume(results['dune_x'][1], initial_vbeach, dclos, dfoot)
    #     initial_bar_vol = results['volbar'][0]
    # else:
    #     # Fix: derive initial vbeach from the initial geometry (results['vbeach'][0]
    #     # is the beach volume *after* the first time step, not the initial state).
    #     initial_vbeach = ((results['dune_x'][0]['shoreline'] - results['dune_x'][0]['ramp_toe']) - bl) / al
    #     initial_dune_vol = calc_dune_vol(results['dune_x'][0], results['dune_z'][0], substrate_x, substrate_z)
    #     initial_beach_vol = calculate_actual_beach_volume(results['dune_x'][0], initial_vbeach, dclos, dfoot)
    #     initial_bar_vol = results['volbar'][0] - results['qbar'][0] * dt  # true initial bar vol (before first update)
    # initial_total = initial_dune_vol + initial_beach_vol + initial_bar_vol
    
    # final_dune_vol = calc_dune_vol(results['dune_x'][-1], results['dune_z'][-1], substrate_x, substrate_z)
    # final_beach_vol = calculate_actual_beach_volume(results['dune_x'][-1], results['vbeach'][-1], dclos, dfoot)
    # final_bar_vol = results['volbar'][-1]
    # final_total = final_dune_vol + final_beach_vol + final_bar_vol
    
    # report_print("INITIAL STATE:")
    # report_print(f"  Dune volume:  {initial_dune_vol:>12.4f} m³/m")
    # report_print(f"  Beach volume: {initial_beach_vol:>12.4f} m³/m")
    # report_print(f"  Bar volume:   {initial_bar_vol:>12.4f} m³/m")
    # report_print(f"  Total volume: {initial_total:>12.4f} m³/m")
    # report_print()

    # report_print("FINAL STATE:")
    # report_print(f"  Dune volume:  {final_dune_vol:>12.4f} m³/m")
    # report_print(f"  Beach volume: {final_beach_vol:>12.4f} m³/m")
    # report_print(f"  Bar volume:   {final_bar_vol:>12.4f} m³/m")
    # report_print(f"  Total volume: {final_total:>12.4f} m³/m")
    # report_print()
    
    # # Calculate changes
    # delta_dune = final_dune_vol - initial_dune_vol
    # delta_beach = final_beach_vol - initial_beach_vol
    # delta_bar = final_bar_vol - initial_bar_vol
    # delta_total = final_total - initial_total
    
    # report_print("VOLUME CHANGES:")
    # report_print(f"  Δ Dune:  {delta_dune:>+12.4f} m³/m ({100*delta_dune/initial_dune_vol:>+7.2f}%)")
    # report_print(f"  Δ Beach: {delta_beach:>+12.4f} m³/m ({100*delta_beach/initial_beach_vol:>+7.2f}%)")
    # report_print(f"  Δ Bar:   {delta_bar:>+12.4f} m³/m ({100*delta_bar/initial_bar_vol:>+7.2f}%)")
    # report_print(f"  Δ Total: {delta_total:>+12.4f} m³/m ({100*delta_total/initial_total:>+7.2f}%)")
    # report_print()
    
    # # Calculate integrated fluxes
    # total_qs_ramp = np.sum(results['qs_ramp'][flux_start:n_steps])
    # total_qs_crest = np.sum(results['qs_crest'][flux_start:n_steps])
    # total_qs_lee = np.sum(results['qs_lee'][flux_start:n_steps])
    # total_aeolian = total_qs_ramp + total_qs_crest + total_qs_lee
    
    # total_qwnet_vol = np.sum([results['qwnet'][i] * dt for i in range(flux_start, n_steps)])
    # total_qbar_vol = np.sum([results['qbar'][i] * dt for i in range(flux_start, n_steps)])
    # total_eroded = np.sum(results['eroded_volume'][flux_start:n_steps])
    
    # flux_label = "over simulation (excluding first step)" if skip_first_step else "over entire simulation"
    # report_print(f"INTEGRATED FLUXES ({flux_label}):")
    # report_print(f"  Aeolian transport from beach to dune: {total_aeolian:>12.4f} m³/m")
    # report_print(f"    - Deposited on ramp:  {total_qs_ramp:>12.4f} m³/m")
    # report_print(f"    - Deposited on crest: {total_qs_crest:>12.4f} m³/m")
    # report_print(f"    - Deposited on lee:   {total_qs_lee:>12.4f} m³/m")
    # report_print(f"  Bar-beach exchange (to bar if +):     {total_qbar_vol:>12.4f} m³/m")
    # report_print(f"  Wave erosion (from dune):             {-total_eroded:>12.4f} m³/m")
    # report_print()
    
    # # Phantom correction: vbeachtot = ramp_toe*dfoot + shoreline*dclos + vbeach.
    # # When the dune toe (ramp_toe) moves by Δrt, vbeachtot changes by
    # # Δrt*(dfoot - db_old) even without any sediment flux.  This geometric
    # # effect is tracked step-by-step and added to the expected beach change so
    # # that the beach error reflects only real flux accounting discrepancies.
    # # When skip_first_step=True, the loop starts at i=flux_start so that the
    # # phantom geometry correction is consistent with the flux integration range.
    # phantom_correction = 0.0
    # vb_prev = initial_vbeach
    # for i in range(flux_start, n_steps):
    #     rt_old = results['dune_x'][i]['ramp_toe']
    #     sl_old = results['dune_x'][i]['shoreline']
    #     rt_new = results['dune_x'][i + 1]['ramp_toe']
    #     beach_width_old = sl_old - rt_old
    #     db_old = vb_prev / beach_width_old if beach_width_old > 0 else 0.0
    #     phantom_correction += (rt_new - rt_old) * (dfoot - db_old)
    #     vb_prev = results['vbeach'][i]
    
    # # Mass balance assessment
    # expected_dune_change = total_aeolian - total_eroded
    # expected_beach_change = -total_aeolian - total_qbar_vol + total_eroded + phantom_correction
    # expected_bar_change = total_qbar_vol
    # # The phantom_correction accounts for the geometric change in vbeachtot when
    # # the dune toe (ramp_toe) moves: the metric vbeachtot = ramp_toe*dfoot +
    # # shoreline*dclos + vbeach shifts even without any real sediment exchange.
    # # This replaces the previous assumption of 0.0.  Since there are no external
    # # sediment sources or sinks, the sum of all three expected compartment
    # # changes must equal phantom_correction (the dune and bar terms cancel out,
    # # leaving only this geometric offset from the beach metric definition).
    # expected_total_change = phantom_correction
    
    # report_print("MASS BALANCE ASSESSMENT:")
    # report_print()
    # report_print("Expected changes (based on fluxes + geometric beach correction):")
    # report_print(f"  Expected Δ Dune:  {expected_dune_change:>+12.4f} m³/m")
    # report_print(f"  Expected Δ Beach: {expected_beach_change:>+12.4f} m³/m")
    # report_print(f"    (incl. phantom correction: {phantom_correction:>+.4f} m³/m)")
    # report_print(f"  Expected Δ Bar:   {expected_bar_change:>+12.4f} m³/m")
    # report_print(f"  Expected Δ Total: {expected_total_change:>+12.4f} m³/m")
    # report_print()
    
    # dune_error = delta_dune - expected_dune_change
    # beach_error = delta_beach - expected_beach_change
    # bar_error = delta_bar - expected_bar_change
    # total_error = delta_total - expected_total_change
    
    # report_print("Mass balance errors:")
    # report_print(f"  Dune error:  {dune_error:>+12.6f} m³/m ({100*abs(dune_error)/initial_dune_vol:.4f}%)")
    # report_print(f"  Beach error: {beach_error:>+12.6f} m³/m ({100*abs(beach_error)/initial_beach_vol:.4f}%)")
    # report_print(f"  Bar error:   {bar_error:>+12.6f} m³/m ({100*abs(bar_error)/initial_bar_vol:.4f}%)")
    # report_print(f"  TOTAL ERROR: {total_error:>+12.6f} m³/m ({100*abs(total_error)/initial_total:.6f}%)")
    # report_print()
    
    # # Time-series conservation analysis (if available)
    # if 'conservation' in results and len(results['conservation'].get('total_error', [])) > 0:
    #     report_print("TIME-SERIES CONSERVATION ANALYSIS:")
    #     report_print()
        
    #     conservation = results['conservation']
        
    #     # Find maximum errors
    #     max_dune_error = np.max(np.abs(conservation['dune_error']))
    #     max_beach_error = np.max(np.abs(conservation['beach_error']))
    #     max_bar_error = np.max(np.abs(conservation['bar_error']))
    #     max_total_error = np.max(np.abs(conservation['total_error']))
        
    #     # Find when maximum errors occurred
    #     idx_max_dune = np.argmax(np.abs(conservation['dune_error']))
    #     idx_max_beach = np.argmax(np.abs(conservation['beach_error']))
    #     idx_max_bar = np.argmax(np.abs(conservation['bar_error']))
    #     idx_max_total = np.argmax(np.abs(conservation['total_error']))
        
    #     # Calculate time for each max error
    #     time_max_dune = idx_max_dune * dt / 3600
    #     time_max_beach = idx_max_beach * dt / 3600
    #     time_max_bar = idx_max_bar * dt / 3600
    #     time_max_total = idx_max_total * dt / 3600
        
    #     report_print("Maximum errors during simulation:")
    #     report_print(f"  Dune:  {max_dune_error:>12.6e} m³/m at t={time_max_dune:>7.2f}h (step {idx_max_dune})")
    #     report_print(f"  Beach: {max_beach_error:>12.6e} m³/m at t={time_max_beach:>7.2f}h (step {idx_max_beach})")
    #     report_print(f"  Bar:   {max_bar_error:>12.6e} m³/m at t={time_max_bar:>7.2f}h (step {idx_max_bar})")
    #     report_print(f"  Total: {max_total_error:>12.6e} m³/m at t={time_max_total:>7.2f}h (step {idx_max_total})")
    #     report_print()
        
    #     # Calculate cumulative absolute errors
    #     cumulative_dune = np.sum(np.abs(conservation['dune_error']))
    #     cumulative_beach = np.sum(np.abs(conservation['beach_error']))
    #     cumulative_bar = np.sum(np.abs(conservation['bar_error']))
    #     cumulative_total = np.sum(np.abs(conservation['total_error']))
        
    #     report_print("Cumulative absolute errors:")
    #     report_print(f"  Dune:  {cumulative_dune:>12.6e} m³/m")
    #     report_print(f"  Beach: {cumulative_beach:>12.6e} m³/m")
    #     report_print(f"  Bar:   {cumulative_bar:>12.6e} m³/m")
    #     report_print(f"  Total: {cumulative_total:>12.6e} m³/m")
    #     report_print()
        
    #     # Identify problematic timesteps (errors > 0.01% of initial volume)
    #     threshold_pct = 0.01
    #     problematic_steps = []
    #     for i in range(len(conservation['total_error'])):
    #         if (conservation['dune_error_pct'][i] > threshold_pct or 
    #             conservation['beach_error_pct'][i] > threshold_pct or 
    #             conservation['bar_error_pct'][i] > threshold_pct):
    #             problematic_steps.append(i)
        
    #     if problematic_steps:
    #         n_show = min(5, len(problematic_steps))
    #         report_print(f"Timesteps with errors > {threshold_pct}% (showing first {n_show} of {len(problematic_steps)}):")
    #         for step in problematic_steps[:n_show]:
    #             time_h = step * dt / 3600
    #             report_print(f"  Step {step:4d} (t={time_h:7.2f}h): "
    #                       f"Dune={conservation['dune_error_pct'][step]:.4f}%, "
    #                       f"Beach={conservation['beach_error_pct'][step]:.4f}%, "
    #                       f"Bar={conservation['bar_error_pct'][step]:.4f}%")
    #     else:
    #         report_print(f"✓ No timesteps with errors > {threshold_pct}% of initial volumes")
        
    #     report_print()
    
    # # Assessment
    # error_threshold = 0.1  # 0.1% threshold
    # total_error_percent = abs(100*total_error/initial_total)
    # conserved = total_error_percent < error_threshold
    
    # if conserved:
    #     report_print("✓ RESULT: Mass is conserved within acceptable limits!")
    #     report_print(f"  (Error < {error_threshold}% threshold)")
    # else:
    #     report_print("⚠ WARNING: Mass balance error exceeds acceptable threshold!")
    #     report_print(f"  (Error > {error_threshold}% threshold)")
    
    # report_print()
    # report_print("=" * 80)
    
    # # Join all report lines into a single string
    # report_text = "\n".join(report_lines)
    
    # return {
    #     'conserved': conserved,
    #     'total_error_percent': total_error_percent,
    #     'total_error': total_error,
    #     'initial_total': initial_total,
    #     'final_total': final_total,
    #     'dune_error': dune_error,
    #     'beach_error': beach_error,
    #     'bar_error': bar_error,
    #     'error_threshold': error_threshold,
    #     'report_text': report_text
    # }
