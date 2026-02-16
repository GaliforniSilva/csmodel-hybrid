"""
main.py
Authors:
    Filipe Galiforni Silva (TU Delft)
    Ad Reniers (TU Delft)
    Sierd de Vries (TU Delft)
    Floris van Rees (Deltares)
    Caroline Hallin (Lund University).

TU Delft, 2025

Description:
    1D Version of the CS Model (Hallin et. al., 2019) adapted to include vegetation and hard structures
    to simulate hybrid dune-dike systems. The model simulates cross-shore profile evolution under the influence
    of waves, water levels, vegetation growth, and the presence of hard structures.

Usage (more info on the README file):
    python main.py path/to/config.yaml

Inputs (within the YAML configuration file):
    - Initial cross-shore profile data
    - Wave and water level data
    - Vegetation parameters
    - Hard structure parameters

Outputs:
    - Updated cross-shore profile over time
    - Updated vegetation coverage over time

Dependencies:
    - numpy
    - scipy
    - matplotlib
    - shapely
    - pandas
    - netCDF4

License:
    Copyright (c) 2025, Delft University of Technology.
    All rights reserved.

References:
HALLIN, C., LARSON, M. & HANSON, H. 2019. Simulating beach and dune evolution at decadal to centennial scale under rising sea levels. PLoS One, 14, e0215651. doi: 10.1371/journal.pone.0215651

"""

############################################################################
# IMPORTS
############################################################################

import sys
import os
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from scipy.interpolate import interp1d
import pandas as pd
import netCDF4 as nc
import yaml
from types import SimpleNamespace


############################################################################
# PATH SETUP
############################################################################

# if len(sys.argv) < 2:
#     raise RuntimeError(
#         "Usage: python main.py path/to/config.yaml"
#     )

#config_path = sys.argv[1]

package_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(package_dir)  # parent folder

if project_root not in sys.path:
    sys.path.insert(0, project_root)

from cs_model.functions import *
import pickle

############################################################################
# INITIALIZATION OF PARAMETERS
############################################################################

def load_config(config_path):
    """
    Load YAML configuration file from a given path.
    """
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    return cfg


def run_model(config_path, check_mass_conservation=False):
    """
    Run the CS model using a YAML configuration file.

    Parameters
    ----------
    config_path : str
        Path to the YAML configuration file
    check_mass_conservation : bool, optional
        If True, perform mass conservation check after simulation and print report.
        Default is False.

    Returns
    -------
    results : dict
        Model output
    """

    # Load configuration
    config = load_config(config_path)


    # Inject each config key as a variable in the global namespace
    # for section, params in config.items():
    #     if isinstance(params, dict):
    #         for k, v in params.items():
    #             globals()[k] = v
    #     else:
    #         globals()[section] = params

    params = {}
    for section, values in config.items():
        if isinstance(values, dict):
            params.update(values)
        else:
            params[section] = values

    p = SimpleNamespace(**params)

    # load boundary conditions
    time_length_dt = int(p.time_length_hours * 3600 / p.dt)

    wave_time_series = load_timeseries(
        p.wave_height, time_length_dt, column="wave_height"
    )

    wave_period_time_series = load_timeseries(
        p.wave_period, time_length_dt, column="wave_period"
    )

    wave_dir_time_series = load_timeseries(
        p.wave_direction, time_length_dt, column="wave_direction"
    )

    swl_time_series = load_timeseries(
        p.swl, time_length_dt, column="swl"
    )

    wind_speed_time_series = load_timeseries(
        p.wind_speed, time_length_dt, column="wind_speed"
    )

    wind_dir_time_series = load_timeseries(
        p.wind_direction, time_length_dt, column="wind_direction"
    )

    month_timeseries = load_timeseries(
        p.month, time_length_dt, column="month", default=7
    )


    # Allocate memory
    qbar = 0#np.zeros(time_length_hours)
    qnour = 0#np.zeros(time_length_hours) #onshore transport from shoreface nourishments
    qslr = 0#np.zeros(time_length_hours)
    dvol = 0#np.zeros(time_length_hours)

    IV = 1 #counter for nourishments
    vnourd = np.zeros(time_length_dt) #nourished volume dune
    vnourb = np.zeros(time_length_dt) #nourished volume beach
    vnourbar = np.zeros(time_length_dt) #nourished volume bar

    volat = 0
    volab = 0
    volad = 0
    volrb = 0

    barnour = np.ones(time_length_dt)*p.barnour_total
    eroded_volume_time_series = np.zeros(time_length_dt)
    overwash_volume_time_series = np.zeros(time_length_dt)
    dune_height_time_series = np.zeros(time_length_dt)
    dune_foot_time_series = np.zeros(time_length_dt)
    total_water_level_time_series = np.zeros(time_length_dt)
    runup_time_series = np.zeros(time_length_dt)
    X_time_series = np.zeros(time_length_dt)
    vbeach_volume_time_series = np.zeros(time_length_dt)

    ##### Define initial profile #####

    # define non-erodible substrate (dike + flat beach base)
    substrate_x = [0, p.yl_hs, p.ylprim_hs, p.ysprim_hs, p.ys_hs, p.yg]
    substrate_z = [p.base_elevation_back, p.base_elevation_back, p.s_hs + p.base_elevation_front, p.s_hs + p.base_elevation_front, p.base_elevation_front, p.base_elevation_front] # use base 0 = mean sea level throughout the profile for the calculations

    # Define dune profile as dictionaries for clarity
    dune_x = {
        'back_ref': 0.0,
        'lee_foot': p.yl,
        'lee_crest': p.ylprim,
        'crest_temp': p.ysprim,
        'seaward_crest': p.ysprim,
        'ramp_dike_temp': p.ysprim, # temporary point to help with calculations
        'ramp_crest': p.ysprim,  # initially same as seaward_crest
        'ramp_toe_temp': p.ys,
        'ramp_toe': p.ys,
        'shoreline': p.yg
    }

    dune_z = {
        'back_ref': p.dfoot,
        'lee_foot': p.dfoot,
        'lee_crest': p.s + p.dfoot,
        'crest_temp': p.s + p.dfoot,
        'seaward_crest': p.s + p.dfoot,
        'ramp_dike_temp': p.sramp + p.dfoot, # temporary point to help with calculations
        'ramp_crest': p.sramp + p.dfoot,
        'ramp_toe_temp': p.dfoot,
        'ramp_toe': p.dfoot,
        'shoreline': 0.0
    }

    # Initialize results dictionary to store time series
    results = {
        'dune_x': [dune_x.copy()],
        'dune_z': [dune_z.copy()],
        'substrate_x': substrate_x,
        'substrate_z': substrate_z,
        'vramp': [],
        'vbeach': [],
        'vdune': [],
        'wave_height': [],
        'swl': [],
        'runup': [],
        'total_water_level': [],
        'eroded_volume': [],
        'overwash_volume': [],
        'dune_height': [dune_z['seaward_crest'] - dune_z['ramp_toe']],
        'dune_foot': [dune_x['ramp_toe']],
        'rhov_df': [],
        'rhov_dt': [],
        'rhov_db': [],
        'rhov_df2': [],
        'rhov_dt2': [],
        'rhov_db2': [],
        'volbar': [],
        'qs_crest': [],
        'qs_lee': [],
        'qs_ramp': [],
        'qwnet': [],
        'qbar': [],
        'wave_period': wave_period_time_series,
        'wave_direction': wave_dir_time_series,
        'wind_speed': wind_speed_time_series,
        'wind_direction': wind_dir_time_series,
        'time': np.arange(0, time_length_dt*p.dt, p.dt)
    }

    # Initialize cumulative time for vegetation update
    cumdt = 0
    vegupdate = 0
    s_veg = p.s # allocate value for the vegetation growth calc
    ys_veg = p.ys # allocate value for the vegetation growth calc
    yl_veg = p.yl # allocate value for the vegetation growth calc
    sramp_veg = p.sramp # allocate value for the vegetation growth calc

    # Initialize varying parameters
    # Calculate initial vbeach from shoreline position (yg) instead of using config value
    # to ensure consistency with the beach volume formula
    vbeach = ((p.yg - p.ys) - p.bl) / p.al
    
    # Initialize ramp volume from config
    # NOTE: vramp must be geometrically consistent with dune parameters
    # Formula: vramp = 0.5 × (ys - ysprim)² × betas
    # Use validate_yaml_volumes.py to verify consistency before running
    vramp = p.vramp
    rhov_df = p.rhov_df
    rhov_dt = p.rhov_dt
    rhov_db = p.rhov_db
    rhov_df2 = p.rhov_df2
    rhov_dt2 = p.rhov_dt2
    rhov_db2 = p.rhov_db2
    volbar = p.volbar
    deltavall = 0

    ############################################################################
    # START SIMULATION
    ############################################################################

    for i in range(len(wave_time_series)):

        month_i = month_timeseries[i]
        ### update counter for vegetation
        cumdt = cumdt + p.dt
        if cumdt >= p.vegdt:
            # Update vegetation parameters here
            cumdt = 0
            vegupdate = 1

        ### hydrodynamics ###
        h0 = wave_time_series[i] / 1.414 #rms wave height
        tp = wave_period_time_series[i]  # Wave period (s)
        swl = swl_time_series[i]  # Sea water level (m)

        # modify wave height with respect to wave direction
        if 90 <= p.shorn <= 270:
            UL = p.shorn + 90
            LL = p.shorn - 90
            if LL <= wave_dir_time_series[i] <= UL:
                teta = wave_dir_time_series[i] - p.shorn
            else:
                teta = 0.0
                h0 = 0.0
                t = 0.0
        elif p.shorn > 270:
            UL = p.shorn + 90 - 360
            LL = p.shorn - 90
            if wave_dir_time_series[i] >= LL:
                teta = wave_dir_time_series[i] - p.shorn
            elif wave_dir_time_series[i] <= UL:
                teta = wave_dir_time_series[i] + 360 - p.shorn
            else:
                teta = 0.0
                h0 = 0.0
                t = 0.0
        elif p.shorn < 90:
            UL = p.shorn + 90
            LL = p.shorn - 90 + 360
            if wave_dir_time_series[i] <= UL:
                teta = wave_dir_time_series[i] - p.shorn
            elif wave_dir_time_series[i] >= LL:
                teta = 360 - wave_dir_time_series[i] + p.shorn
            else:
                teta = 0.0
                h0 = 0.0
                tp = 0.0

        # Do not allow the waves to become too small
        if tp <= 0.0:
            tp = 1.0
        if h0 <= 0.0:
            h0 = 0.1

        # calculate approximate runup
        
        LO = 1.5613 * (0.95 * tp) ** 2
        slope_foreshore = dune_z['ramp_toe']/(dune_x['shoreline']- dune_x['ramp_toe'])
        r = (0.158/0.03)*slope_foreshore*np.sqrt(h0 * LO) #we divide by an actual slope and recalculate based on the slope we have so to vary potential dissipation with X
        X = 2 * vbeach / p.dfoot * (1 - swl / p.dfoot)

        if X < 0:
            X = 0 # surge above the dunefoot, hence no bed for runup friction, only wave dissipation.
            
        rprim = (r) * np.exp(-2 * p.cf * X) + (p.dfoot) * (1 - np.exp(-2 * p.cf * X))
            
        if rprim > p.dfoot:
            r = rprim

        runup_time_series[i] = r
        X_time_series[i] = X
        
        aover = 3.0
        alfa = 0.0
        # calculate potential erosion

        # total frontal area per alongshore cm (cm^2/cm)
        Av_cm2_per_cm_spec1 = max(0.0, p.frontal_cover_scale_spec1*rhov_df + p.frontal_cover_scale_spec1*rhov_dt + p.frontal_cover_scale_spec1*rhov_db)
        Av_cm2_per_cm_spec2 = max(0.0, p.frontal_cover_scale_spec2*rhov_df2 + p.frontal_cover_scale_spec2*rhov_dt2 + p.frontal_cover_scale_spec2*rhov_db2)

        Av_total = Av_cm2_per_cm_spec1 + Av_cm2_per_cm_spec2
        # Figlus-like slope: 2.1% reduction per (cm^2/cm) -> 0.021 fraction per (cm^2/cm)
        # multiplier to apply to qd (1 = no effect, 0.72 = 28% reduction)
        vegred = (1.0 - min(0.021 * Av_total, 0.28))

        if r + swl > p.dfoot:
            if r + swl <= dune_z['seaward_crest']:
                qd = vegred * 4 * p.cimpact * (r + swl - p.dfoot) ** 2 / tp
            else:
                print('OVERWASH')
                qd = vegred * 4 * p.cimpact * (((r + swl - p.dfoot) * (dune_z['seaward_crest']-p.dfoot)) / tp)
                alfa = (((r + swl - p.dfoot) / (dune_z['seaward_crest']-p.dfoot)) - 1) / aover
        else:
            qd = 0.0
        
        qs = qd / (1 + alfa)
        ql = qd * alfa / (1 + alfa)

        total_water_level = r + swl 

        total_water_level_time_series[i] = total_water_level

        #### aeolian transport at the beach ####
            ## calculate wind direction to the beach normal
        if 80 <= p.shorn <= 280:
            UL = p.shorn + 80
            LL = p.shorn - 80
            if LL <= wind_dir_time_series[i] <= UL:
                tetaw = wind_dir_time_series[i] - p.shorn
                wspeed = wind_speed_time_series[i]
            else:
                tetaw = 0.0
                wspeed = 0.0
        elif p.shorn > 280:
            UL = p.shorn + 80 - 360
            LL = p.shorn - 80
            if wind_dir_time_series[i] >= LL:
                tetaw = wind_dir_time_series[i] - p.shorn
                wspeed = wind_speed_time_series[i]
            elif wind_dir_time_series[i] <= UL:
                tetaw = wind_dir_time_series[i] + 360 - p.shorn
                wspeed = wind_speed_time_series[i]
            else:
                tetaw = 0.0
                wspeed = 0.0
        elif p.shorn < 80:
            UL = p.shorn + 80
            LL = p.shorn - 80 + 360
            if wind_dir_time_series[i] <= UL:
                tetaw = wind_dir_time_series[i] - p.shorn
                wspeed = wind_speed_time_series[i]
            elif wind_dir_time_series[i] >= LL:
                tetaw = 360 - wind_dir_time_series[i] + p.shorn
                wspeed = wind_speed_time_series[i]
            else:
                tetaw = 0.0
                wspeed = 0.0
        
        wdir = tetaw # final adjusted wind direction

        ## Calculate wind shear and critical shear at the beach
        wcrit = p.awind * np.sqrt((p.rhos - p.rhoa) / p.rhoa * p.g * p.d50)
        #z0 = 0.081 * np.log10(d50 * 10 ** 3 / 0.18) # the one used in Hallin's model
        z0 = p.d50/30 # the one used in Bagnold

        #wshear = wspeed * max(kappa / np.log(z / z0), 0.053)
        wshear = wspeed * p.kappa / np.log(p.z / z0)  # shear velocity

        if wshear > wcrit:
            mwe = p.kw * np.sqrt(p.d50 / p.d50ref) * p.rhoa * wshear ** 2 / p.g * (wshear - wcrit)
        else:
            mwe = 0
        qwe = mwe / (p.rhos * (1 - p.p))

        #yr = ys + (1 - (r + swl) / dfoot) * (yg - ys)
        yr = dune_x['ramp_toe'] + (1 - (r + swl) / p.dfoot) * (dune_x['shoreline'] - dune_x['ramp_toe'])
        #bdry = yr - ys
        bdry = yr - dune_x['ramp_toe']

        if r + swl >= p.dfoot or yr < dune_x['ramp_toe']:
            bdry = 0
        f = bdry / abs(np.cos(np.radians(wdir)))
        qw = qwe * (1 - np.exp(-p.cwrate * f))
        qwnet = qw * np.cos(np.radians(wdir))

        # Introduce dynamic reduction based on Duran and Moore (2013) Ad
        #qwnet = qwnet * max(1.0 - 4.2 * (s-p.dfoot) / (yg - ys), 0.0)
        qwnet = qwnet * max(1.0 - 4.2 * (dune_z['seaward_crest']-p.dfoot) / (dune_x['shoreline'] - dune_x['ramp_toe']), 0.0)

        ######## AEOLIAN DEPOSITION AT DUNE ########

        qs_ramp, qs_crest, qs_lee = deposition_dune(dune_x, dune_z, substrate_x, substrate_z, qwnet, p.dt, tetaw, wspeed, p.rhos, p.rhoa, p.d50, p.kappa, p.z, p.vegcount, p.awind, rhov_df=rhov_df, rhov_dt=rhov_dt, rhov_db=rhov_db, 
            rhov_df2=rhov_df2 if p.vegcount == 1 else 0, 
            rhov_dt2=rhov_dt2 if p.vegcount == 1 else 0, 
            rhov_db2=rhov_db2 if p.vegcount == 1 else 0)

        ######## NOURISHMENT AND SUBTIDAL EVOLUTION ########
        
        wfall = FALVEL(p.temp, p.d50) # calculate fall velocity
        
        # subtidal evolution
        qbar, qnour, qslr, dvol = subaq(p.dt, volbar, h0, tp, wfall, p.slr, barnour[i], p.cb, p.bactive)

        # specify nourished volumes
        if p.vanr > 0:
            vnourd[i], vnourb[i], vnourbar[i], IV = nour(i, IV, volat, volad, volab, volrb, p.vanr)
            
        # Calculate bar evolution
        volbar, barnour[i] = bar(volbar, qbar, p.dt, vnourbar[i], barnour[i], qnour)
        
        # Calculate beach volume BEFORE morpho_update (using old ramp_toe)
        vbeach_old = ((dune_x['shoreline'] - dune_x['ramp_toe'])-p.bl)/p.al
        db_old = vbeach_old / (dune_x['shoreline'] - dune_x['ramp_toe'])
        vbeachtot_old = dune_x['ramp_toe']  * (p.dclos + p.dfoot) + (dune_x['shoreline'] - dune_x['ramp_toe']) * (p.dclos + db_old)
        
        if i==51568:
            print(i)
        ####### MORPHOLOGICAL UPDATE OF THE DUNE PROFILE #########
        dune_x, dune_z, deltavramp, erod_overwash, deltavramp_ero, deltavall = morpho_update(dune_x, dune_z, substrate_x, substrate_z, qs_ramp, qs_crest, qs_lee, qs, ql, r, swl, vnourd[i], p.dt, p.betas, p.betas_hs, p.betal, p.betal_hs, vramp, deltavall)

        ######## VEGETATION GROWTH UPDATE #########
            # update vegetation parameters
        if p.vegcount == 1 and vegupdate == 1:
            rhov_df, rhov_dt, rhov_db, rhov_df2, rhov_dt2, rhov_db2 = vegNield_seasonal(p.vegdt,s_veg,dune_z['seaward_crest'], dune_x['lee_foot'], yl_veg, dune_x['ramp_toe'], ys_veg, dune_z['ramp_crest'], sramp_veg, p.betas, p.betal, rhov_df, rhov_dt, rhov_db, rhov_df2, rhov_dt2, rhov_db2, month_i)
            vegupdate = 0
            s_veg = dune_z['seaward_crest']
            ys_veg = dune_x['ramp_toe']
            yl_veg = dune_x['lee_foot']
            sramp_veg = dune_z['ramp_crest']
        elif p.vegcount == 0 and vegupdate == 1:
            rhov_df, rhov_dt, rhov_db = veg(ys, ysnew, yl, ylnew, sramp, srampnew, betas, betal, s, snew, vegdt,rhov_df,rhov_dt,rhov_db)
            vegupdate = 0

        
        # update values for next iteration
        vramp = partial_volume(dune_x, dune_z, substrate_x, substrate_z, dune_x['ramp_crest'], dune_x['ramp_toe'])
        vdune = partial_volume(dune_x, dune_z, substrate_x, substrate_z, dune_x['lee_crest'], dune_x['ramp_toe'])
        dune_height_time_series[i] = dune_z['seaward_crest'] - dune_z['ramp_toe']
        dune_foot_time_series[i] = dune_x['ramp_toe']
        total_water_level_time_series[i] = total_water_level
        eroded_volume_time_series[i] = deltavramp #erod_ramp
        overwash_volume_time_series[i] = erod_overwash
        
        # Calculate new beach volume including eroded sediment from dune
        vbeachtotnew = vbeachtot_old + (-qwnet - qslr + qnour - qbar) * p.dt + vnourb[i] + deltavramp_ero
        vbeachtotmin = dune_x['ramp_toe'] * (p.dclos + p.dfoot) + 0.5 * p.dfoot**2 / p.betaf + p.dfoot / p.betaf * p.dclos

        # Calculate new shoreline coordinate
        dune_x['shoreline'] = dune_x['ramp_toe'] + (p.al * (vbeachtotnew - dune_x['ramp_toe'] * (p.dfoot + p.dclos)) + p.bl) / (1 + p.al * p.dclos)
        
        # Calculate final beach volume for output
        vbeach = ((dune_x['shoreline'] - dune_x['ramp_toe'])-p.bl)/p.al
        vbeach_volume_time_series[i] = vbeach
        
        ####### STORE OUTPUTS #######
        results['dune_x'].append(dune_x.copy())
        results['dune_z'].append(dune_z.copy())
        results['vramp'].append(vramp)
        results['vbeach'].append(vbeach)
        results['vdune'].append(vdune)
        results['wave_height'].append(wave_time_series[i])
        results['swl'].append(swl)
        results['runup'].append(r)
        results['total_water_level'].append(total_water_level)
        results['eroded_volume'].append(deltavramp_ero)
        results['overwash_volume'].append(erod_overwash)
        results['dune_height'].append(dune_z['seaward_crest'] - dune_z['ramp_toe'])
        results['dune_foot'].append(dune_x['ramp_toe'])
        results['rhov_db'].append(rhov_db)
        results['rhov_dt'].append(rhov_dt)
        results['rhov_df'].append(rhov_df)
        results['rhov_db2'].append(rhov_db2)
        results['rhov_dt2'].append(rhov_dt2)
        results['rhov_df2'].append(rhov_df2)
        results['volbar'].append(volbar)
        results['qs_crest'].append(qs_crest)
        results['qs_lee'].append(qs_lee)
        results['qs_ramp'].append(qs_ramp)
        results['qwnet'].append(qwnet)
        results['qbar'].append(qbar)

        ## print progress
        if i % 1000 == 0:
            print(f"Progress: {100*i/len(wave_time_series):.1f}% ({i} of {len(wave_time_series)} time steps)")

    ############################################################################
    # SAVE RESULTS
    ############################################################################

    output_filename = config.get('output', {}).get('filename', 'cs_model_results.pkl')
    with open(output_filename, 'wb') as f:
        pickle.dump(results, f)

    print(f"Simulation completed. Results saved to {output_filename}.")
    
    ############################################################################
    # MASS CONSERVATION CHECK (if requested)
    ############################################################################
    
    if check_mass_conservation:
        print()  # Empty line for better formatting
        from cs_model.mass_conservation import check_mass_conservation as check_mass
        mass_check_result = check_mass(results, config, verbose=True)
        
        # Save the full mass conservation report to file
        if mass_check_result:
            mass_report_filename = output_filename.replace('.pkl', '_mass_conservation.txt')
            with open(mass_report_filename, 'w') as f:
                f.write("Detailed Mass Conservation Report\n")
                f.write("=" * 80 + "\n")
                f.write(f"Output file: {output_filename}\n")
                f.write("=" * 80 + "\n")
                f.write("\n")
                # Write the full detailed report
                f.write(mass_check_result['report_text'])
            print(f"Mass conservation report saved to {mass_report_filename}")

            # Print final summary to terminal for visibility
            print()
            print("=" * 80)
            print("MASS CONSERVATION SUMMARY")
            print("=" * 80)
            print(f"Status: {'✓ PASSED' if mass_check_result['conserved'] else '⚠ FAILED'}")
            print(f"Total Volume Error: {mass_check_result['total_error']:+.6f} m³/m ({mass_check_result['total_error_percent']:.6f}%)")
            print(f"Initial Total Volume: {mass_check_result['initial_total']:.4f} m³/m")
            print(f"Final Total Volume: {mass_check_result['final_total']:.4f} m³/m")
            print(f"Error Threshold: {mass_check_result['error_threshold']:.1f}%")
            if mass_check_result['conserved']:
                print("\n✓ Mass is conserved within acceptable limits!")
            else:
                print("\n⚠ WARNING: Mass balance error exceeds acceptable threshold!")
            print("=" * 80)
            print()

    return results

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Usage: python main.py path/to/config.yaml")

    run_model(sys.argv[1])
