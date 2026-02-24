'''
functions.py
Authors:
    Filipe Galiforni Silva (TU Delft)
    Ad Reniers (TU Delft)
    Sierd de Vries (TU Delft)
    Floris van Rees (Deltares)
    Caroline Hallin (Lund University).

TU Delft, 2025

Description:
    This module contains helper functions used in the main.py file.

License:
    Copyright (c) 2025, Delft University of Technology.
    All rights reserved.
'''

############################################################################
# IMPORTS
############################################################################

import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from scipy.interpolate import interp1d
import pandas as pd
import netCDF4 as nc
import yaml

############################################################################
# FUNCTIONS
############################################################################

### General helper functions
def load_timeseries(
    source,
    length,
    column,
    default=None,
    sep="\t",
    time_col="datetime"
):
    """
    Load a time series from file or create a constant array.

    Parameters
    ----------
    source : str or float
        Filename or constant value.
    length : int
        Length of the output time series.
    column : str
        Column name in file to read.
    default : float, optional
        Default constant value if source is None.
    sep : str
        File separator.
    time_col : str
        Name of datetime column.

    Returns
    -------
    np.ndarray
    """
    if isinstance(source, str):
        df = pd.read_csv(source, sep=sep, parse_dates=[time_col])
        df.set_index(time_col, inplace=True)
        return df[column].values

    if source is None:
        if default is None:
            raise ValueError("source is None and no default provided")
        return np.full(length, default)

    return np.full(length, source)


def trapezoid_positive(y, x):
    """
    Integrate using trapezoid, counting only positive values.
    Negative values are clipped to zero.
    """
    y_clip = np.maximum(y, 0)
    # Use trapezoid (numpy 2.x) or trapz (numpy 1.x) for compatibility
    try:
        return np.trapezoid(y_clip, x)
    except AttributeError:
        return np.trapz(y_clip, x)

def load_config(path="config.yaml"):
    with open(path, "r") as f:
        return yaml.safe_load(f)

### Functions for volume calculations
def ramp_volume(sramp_potential, sramp, ysprim, ys, betas, dfoot, substrate_x, substrate_z):
    '''
    Calculate volume of the ramp part of the dune.
    '''
    # toe position consistent with slope
    yspotential = ys + ((sramp_potential-sramp) / betas)

    # ramp line
    dune_x_line = np.array([ysprim, yspotential])
    dune_z_line = np.array([sramp_potential + dfoot, dfoot])
    
    n_per_segment = 100  # how many points per segment
    x_dense = np.concatenate([np.linspace(dune_x_line[i], dune_x_line[i+1], n_per_segment, endpoint=False) for i in range(len(dune_x_line)-1)] + [[dune_x_line[-1]]])
    z_dense = np.concatenate([np.linspace(dune_z_line[i], dune_z_line[i+1], n_per_segment, endpoint=False) for i in range(len(dune_z_line)-1)] + [[dune_z_line[-1]]])

    # substrate under ramp
    z_sub_line = np.interp(x_dense, substrate_x, substrate_z)

    # trapezoidal integration
    return trapezoid_positive(z_dense - z_sub_line, x_dense)

def lee_volume(ylprim_potential, dune_x, dune_z, betal, substrate_x, substrate_z):
    '''
    Calculate volume of the lee part of the dune.
    '''
    # toe position consistent with slope
    #ylpotential = ylprim_potential - ((s_new_c) / betal) 
    ylpotential = ylprim_potential - ((dune_z['lee_crest'] -dune_z['lee_foot']) / betal) 

    # lee line
    #dune_x_line = np.array([ylpotential, ylprim_potential, ylprim_new_c]) #to calculate, keep the initial ylprim always otherwise u will have the same lee always!
    dune_x_line = np.array([ylpotential, ylprim_potential, dune_x['lee_crest']]) #to calculate, keep the initial ylprim always otherwise u will have the same lee always!
    dune_z_line = np.array([dune_z['lee_foot'], dune_z['lee_crest'], dune_z['lee_crest']]) # new dune height

    n_per_segment = 100  # how many points per segment
    x_dense = np.concatenate([np.linspace(dune_x_line[i], dune_x_line[i+1], n_per_segment, endpoint=False) for i in range(len(dune_x_line)-1)] + [[dune_x_line[-1]]])
    z_dense = np.concatenate([np.linspace(dune_z_line[i], dune_z_line[i+1], n_per_segment, endpoint=False) for i in range(len(dune_z_line)-1)] + [[dune_z_line[-1]]])
   
    # substrate under lee side
    z_sub_line = np.interp(x_dense, substrate_x, substrate_z)

    # trapezoidal integration
    return trapezoid_positive(z_dense - z_sub_line, x_dense)

def total_volume(dune_x, dune_z, substrate_x, substrate_z):
    '''
    Compute the total dune volume between dune and substrate.
    Interpolates the dune curve so integration is more accurate.
    '''
    # Remove any None values from dune_x and dune_z (and keep them in sync)
    dune_x, dune_z = zip(*[(x, z) for x, z in zip(dune_x, dune_z) if x is not None and z is not None])
    dune_x = list(dune_x)
    dune_z = list(dune_z)

    n_per_segment = 100  # how many points per segment
    x_dense = np.concatenate([np.linspace(dune_x[i], dune_x[i+1], n_per_segment, endpoint=False) for i in range(len(dune_x)-1)] + [[dune_x[-1]]])
    z_dense = np.concatenate([np.linspace(dune_z[i], dune_z[i+1], n_per_segment, endpoint=False) for i in range(len(dune_z)-1)] + [[dune_z[-1]]])

    # interpolate substrate under dune
    z_sub = np.interp(x_dense, substrate_x, substrate_z)

    # trapezoidal integration
    return trapezoid_positive(z_dense - z_sub, x_dense)

def partial_volume(dune_x, dune_z, substrate_x, substrate_z,
                   x_min=None, x_max=None):

    # Remove None values
    dune_x = {k: v for k, v in dune_x.items() if v is not None}
    dune_z = {k: v for k, v in dune_z.items() if v is not None}

    # --- Convert dicts to ordered arrays ---
    order = [
        'back_ref', 'lee_foot', 'lee_crest', 'crest_temp',
        'seaward_crest', 'ramp_dike_temp', 'ramp_crest',
        'ramp_toe_temp', 'ramp_toe', 'shoreline'
    ]

    x = np.array([dune_x[k] for k in order if k in dune_x])
    z = np.array([dune_z[k] for k in order if k in dune_z])

    # Densify dune curve
    n_per_segment = 100

    x_dense = np.concatenate([
        np.linspace(x[i], x[i + 1], n_per_segment, endpoint=False)
        for i in range(len(x) - 1)
    ] + [[x[-1]]])

    z_dense = np.concatenate([
        np.linspace(z[i], z[i + 1], n_per_segment, endpoint=False)
        for i in range(len(z) - 1)
    ] + [[z[-1]]])

    # Interpolate substrate
    z_sub = np.interp(x_dense, substrate_x, substrate_z)

    # Restrict to x limits
    if x_min is not None or x_max is not None:
        x_min = x_dense[0] if x_min is None else max(x_min, x_dense[0])
        x_max = x_dense[-1] if x_max is None else min(x_max, x_dense[-1])

        z_min = np.interp(x_min, x_dense, z_dense)
        z_max = np.interp(x_max, x_dense, z_dense)
        z_sub_min = np.interp(x_min, x_dense, z_sub)
        z_sub_max = np.interp(x_max, x_dense, z_sub)

        mask = (x_dense >= x_min) & (x_dense <= x_max)
        x_clip = x_dense[mask]
        z_clip = z_dense[mask]
        z_sub_clip = z_sub[mask]

        if len(x_clip) == 0 or x_clip[0] > x_min:
            x_clip = np.insert(x_clip, 0, x_min)
            z_clip = np.insert(z_clip, 0, z_min)
            z_sub_clip = np.insert(z_sub_clip, 0, z_sub_min)

        if x_clip[-1] < x_max:
            x_clip = np.append(x_clip, x_max)
            z_clip = np.append(z_clip, z_max)
            z_sub_clip = np.append(z_sub_clip, z_sub_max)
    else:
        x_clip, z_clip, z_sub_clip = x_dense, z_dense, z_sub

    return trapezoid_positive(z_clip - z_sub_clip, x_clip)


### Vegetation-related functions
def vegNield(dt,sold,snew, ylnew, ylold, ysnew, ysold, srampnew, srampold, betas, betal, rhov_df_spec1, rhov_dt_spec1, rhov_db_spec1, rhov_df_spec2, rhov_dt_spec2, rhov_db_spec2):
    '''
    Vegetation growth function adapted from Nield and Baas, 2008 and DUBEVEG model (Keijsers et. al., 2016), for two species
    '''

    # define vertices values for both species
    sp1_peak = 0.2 # growth peak for species 1
    sp2_peak = 0.05 # growth peak for species 2

    # Convert time step to years
    dtYear = dt/(3600*24*365) # convert time step to years

    # species 1 - pioneer vegetation (e.g. Ammophila arenaria)
    
    # sedimentation
    x1_spec1 = -0.6*dtYear
    x2_spec1 = 0.2*dtYear
    x3_spec1 = 0.6*dtYear
    x4_spec1 = 2.0*dtYear
    x5_spec1 = 2.2*dtYear
    
    #growth
    y1_spec1 = -0.5*dtYear
    y2_spec1 = 0.0*dtYear
    y3_spec1 = sp1_peak*dtYear 
    y4_spec1 = 0.0*dtYear
    y5_spec1 = -0.5*dtYear

    # species 2 - mature vegetation (e.g. Elymus farctus or Buckthrorn)
    x1_spec2 = -1.4*dtYear
    x2_spec2 = -0.65*dtYear
    x3_spec2 = 0.0*dtYear
    x4_spec2 = 0.2*dtYear
    x5_spec2 = 2.8*dtYear

    y1_spec2 = -0.5*dtYear
    y2_spec2 = 0.0*dtYear
    y3_spec2 = sp2_peak*dtYear
    y4_spec2 = 0.0*dtYear
    y5_spec2 = -0.5*dtYear

    #calculate the sedimentation map
    dh_df = 0.0
    if (abs(srampnew-srampold)>0.0):
        dh_df = srampnew-srampold
    elif (abs(ysnew-ysold)>0.0): #movement of the dunetoe
        dh_df=(ysnew-ysold)*betas
            
    dh_dt = (snew-sold)
    dh_db = -(ylnew-ylold)*betal

    # calculate slopes between vertices
    s12_spec1=(y2_spec1-y1_spec1)/(x2_spec1-x1_spec1)
    s23_spec1=(y3_spec1-y2_spec1)/(x3_spec1-x2_spec1)
    s34_spec1=(y4_spec1-y3_spec1)/(x4_spec1-x3_spec1)
    s45_spec1=(y5_spec1-y4_spec1)/(x5_spec1-x4_spec1)

    s12_spec2=(y2_spec2-y1_spec2)/(x2_spec2-x1_spec2)
    s23_spec2=(y3_spec2-y2_spec2)/(x3_spec2-x2_spec2)
    s34_spec2=(y4_spec2-y3_spec2)/(x4_spec2-x3_spec2)
    s45_spec2=(y5_spec2-y4_spec2)/(x5_spec2-x4_spec2)

    # calculate the new vegetation height for dune front for species 1
    if (dh_df<x1_spec1):
        rhov_dfnew_spec1 = 0.0001
    elif (dh_df>=x1_spec1 and dh_df<x2_spec1):
        rhov_dfnew_spec1 = rhov_df_spec1 + ((dh_df-x1_spec1)*s12_spec1 + y1_spec1)
    elif (dh_df>=x2_spec1 and dh_df<x3_spec1):
        rhov_dfnew_spec1 = rhov_df_spec1 + ((dh_df-x2_spec1)*s23_spec1 + y2_spec1)
    elif (dh_df>=x3_spec1 and dh_df<x4_spec1):
        rhov_dfnew_spec1 = rhov_df_spec1 + ((dh_df-x3_spec1)*s34_spec1 + y3_spec1)
    elif (dh_df>=x4_spec1 and dh_df<x5_spec1):
        rhov_dfnew_spec1 = rhov_df_spec1 + ((dh_df-x4_spec1)*s45_spec1 + y4_spec1)
    elif (dh_df>=x5_spec1):
        rhov_dfnew_spec1 = rhov_df_spec1 + y5_spec1

    # calculate the new vegetation height for dune top for species 1
    if (dh_dt<x1_spec1):
        rhov_dtnew_spec1 = 0.0
    elif (dh_dt>=x1_spec1 and dh_dt<x2_spec1):
        rhov_dtnew_spec1 = rhov_dt_spec1 + ((dh_dt-x1_spec1)*s12_spec1 + y1_spec1)
    elif (dh_dt>=x2_spec1 and dh_dt<x3_spec1):
        rhov_dtnew_spec1 = rhov_dt_spec1 + ((dh_dt-x2_spec1)*s23_spec1 + y2_spec1)
    elif (dh_dt>=x3_spec1 and dh_dt<x4_spec1):
        rhov_dtnew_spec1 = rhov_dt_spec1 + ((dh_dt-x3_spec1)*s34_spec1 + y3_spec1)
    elif (dh_dt>=x4_spec1 and dh_dt<x5_spec1):
        rhov_dtnew_spec1 = rhov_dt_spec1 + ((dh_dt-x4_spec1)*s45_spec1 + y4_spec1)
    elif (dh_dt>=x5_spec1):
        rhov_dtnew_spec1 = rhov_dt_spec1 + y5_spec1

    # calculate the new vegetation height for dune back for species 1
    if (dh_db<x1_spec1):
        rhov_dbnew_spec1 = 0.0
    elif (dh_db>=x1_spec1 and dh_db<x2_spec1):
        rhov_dbnew_spec1 = rhov_db_spec1 + ((dh_db-x1_spec1)*s12_spec1 + y1_spec1)
    elif (dh_db>=x2_spec1 and dh_db<x3_spec1):
        rhov_dbnew_spec1 = rhov_db_spec1 + ((dh_db-x2_spec1)*s23_spec1 + y2_spec1)
    elif (dh_db>=x3_spec1 and dh_db<x4_spec1):
        rhov_dbnew_spec1 = rhov_db_spec1 + ((dh_db-x3_spec1)*s34_spec1 + y3_spec1)
    elif (dh_db>=x4_spec1 and dh_db<x5_spec1):
        rhov_dbnew_spec1 = rhov_db_spec1 + ((dh_db-x4_spec1)*s45_spec1 + y4_spec1)
    elif (dh_db>=x5_spec1):
        rhov_dbnew_spec1 = rhov_db_spec1 + y5_spec1

    # calculate the new vegetation height for dune front for species 2 only if there are species 1 or 2 already there

    if rhov_dfnew_spec1 > 0.5 or rhov_df_spec2 > 0.0:
        if (dh_df<x1_spec2):
            rhov_dfnew_spec2 = 0.0
        elif (dh_df>=x1_spec2 and dh_df<x2_spec2):
            rhov_dfnew_spec2 = rhov_df_spec2 + ((dh_df-x1_spec2)*s12_spec2 + y1_spec2)
        elif (dh_df>=x2_spec2 and dh_df<x3_spec2):
            rhov_dfnew_spec2 = rhov_df_spec2 + ((dh_df-x2_spec2)*s23_spec2 + y2_spec2)
        elif (dh_df>=x3_spec2 and dh_df<x4_spec2):
            rhov_dfnew_spec2 = rhov_df_spec2 + ((dh_df-x3_spec2)*s34_spec2 + y3_spec2)
        elif (dh_df>=x4_spec2 and dh_df<x5_spec2):
            rhov_dfnew_spec2 = rhov_df_spec2 + ((dh_df-x4_spec2)*s45_spec2 + y4_spec2)
        elif (dh_df>=x5_spec2):
            rhov_dfnew_spec2 = rhov_df_spec2 + y5_spec2
    else:
        rhov_dfnew_spec2 = 0.0

    # calculate the new vegetation height for dune top for species 2
    if rhov_dtnew_spec1 > 0.5 or rhov_dt_spec2 > 0.0:
        if (dh_dt<x1_spec2):
            rhov_dtnew_spec2 = 0.0
        elif (dh_dt>=x1_spec2 and dh_dt<x2_spec2):
            rhov_dtnew_spec2 = rhov_dt_spec2 + ((dh_dt-x1_spec2)*s12_spec2 + y1_spec2)
        elif (dh_dt>=x2_spec2 and dh_dt<x3_spec2):
            rhov_dtnew_spec2 = rhov_dt_spec2 + ((dh_dt-x2_spec2)*s23_spec2 + y2_spec2)
        elif (dh_dt>=x3_spec2 and dh_dt<x4_spec2):
            rhov_dtnew_spec2 = rhov_dt_spec2 + ((dh_dt-x3_spec2)*s34_spec2 + y3_spec2)
        elif (dh_dt>=x4_spec2 and dh_dt<x5_spec2):
            rhov_dtnew_spec2 = rhov_dt_spec2 + ((dh_dt-x4_spec2)*s45_spec2 + y4_spec2)
        elif (dh_dt>=x5_spec2):
            rhov_dtnew_spec2 = rhov_dt_spec2 + y5_spec2
    else:
        rhov_dtnew_spec2 = 0.0

    # calculate the new vegetation height for dune back for species 2

    if rhov_dbnew_spec1 > 0.5 or rhov_db_spec2 > 0.0:
        if (dh_db<x1_spec2):
            rhov_dbnew_spec2 = 0.0
        elif (dh_db>=x1_spec2 and dh_db<x2_spec2):
            rhov_dbnew_spec2 = rhov_db_spec2 + ((dh_db-x1_spec2)*s12_spec2 + y1_spec2)
        elif (dh_db>=x2_spec2 and dh_db<x3_spec2):
            rhov_dbnew_spec2 = rhov_db_spec2 + ((dh_db-x2_spec2)*s23_spec2 + y2_spec2)
        elif (dh_db>=x3_spec2 and dh_db<x4_spec2):
            rhov_dbnew_spec2 = rhov_db_spec2 + ((dh_db-x3_spec2)*s34_spec2 + y3_spec2)
        elif (dh_db>=x4_spec2 and dh_db<x5_spec2):
            rhov_dbnew_spec2 = rhov_db_spec2 + ((dh_db-x4_spec2)*s45_spec2 + y4_spec2)
        elif (dh_db>=x5_spec2):
            rhov_dbnew_spec2 = rhov_db_spec2 + y5_spec2
    else:   
        rhov_dbnew_spec2 = 0.0

    # cap the vegetation height to a minimum value
    if (rhov_dfnew_spec1 < 0):
        rhov_dfnew_spec1 = 0.0
    if (rhov_dtnew_spec1 < 0):
        rhov_dtnew_spec1 = 0.0
    if (rhov_dbnew_spec1 < 0):
        rhov_dbnew_spec1 = 0.0

    if (rhov_dfnew_spec2 < 0):
        rhov_dfnew_spec2 = 0.0
    if (rhov_dtnew_spec2 < 0):
        rhov_dtnew_spec2 = 0.0
    if (rhov_dbnew_spec2 < 0):
        rhov_dbnew_spec2 = 0.0

    if (rhov_dfnew_spec1 > 1):
        rhov_dfnew_spec1 = 1.0
    if (rhov_dtnew_spec1 > 1):
        rhov_dtnew_spec1 = 1.0
    if (rhov_dbnew_spec1 > 1):
        rhov_dbnew_spec1 = 1.0

    if (rhov_dfnew_spec2 > 1):
        rhov_dfnew_spec2 = 1.0
    if (rhov_dtnew_spec2 > 1):
        rhov_dtnew_spec2 = 1.0
    if (rhov_dbnew_spec2 > 1):
        rhov_dbnew_spec2 = 1.0

    # calculate the balance between the two species
    if (rhov_dfnew_spec1 + rhov_dfnew_spec2 > 1.0): # if the sum of the two species is greater than 1.0, then the difference is taken from the pioneer species
        rhov_dfnew_spec1 = 1.0 - rhov_dfnew_spec2

    if (rhov_dtnew_spec1 + rhov_dtnew_spec2 > 1.0): # if the sum of the two species is greater than 1.0, then the difference is taken from the pioneer species
        rhov_dtnew_spec1 = 1.0 - rhov_dtnew_spec2

    if (rhov_dbnew_spec1 + rhov_dbnew_spec2 > 1.0): # if the sum of the two species is greater than 1.0, then the difference is taken from the pioneer species
        rhov_dbnew_spec1 = 1.0 - rhov_dbnew_spec2

    return rhov_dfnew_spec1, rhov_dtnew_spec1, rhov_dbnew_spec1, rhov_dfnew_spec2, rhov_dtnew_spec2, rhov_dbnew_spec2

def seasonal_factor(month):
    """Returns seasonal growth multiplier (0–1) for pioneer and mature vegetation."""
    # Example: growth in spring/summer, zero in winter
    month_factors = [0.05, 0.05, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.05]

    return month_factors[month-1]

def vegNield_seasonal(dt, sold, snew, ylnew, ylold, ysnew, ysold,
                      srampnew, srampold, betas, betal,
                      rhov_df_spec1, rhov_dt_spec1, rhov_db_spec1,
                      rhov_df_spec2, rhov_dt_spec2, rhov_db_spec2,
                      month=None):
    
    '''
    Vegetation growth function adapted from Nield and Baas, 2008 and DUBEVEG model (Keijsers et. al., 2016), for two species
    with seasonal growth modulation.
    '''
    
    # define vertices values for both species
    sp1_peak = 0.2  # pioneer species peak growth
    sp2_peak = 0.05 # mature species peak growth

    dtYear = dt/(3600*24*365)  # timestep in years
    dtMonth = dt/(3600*24*30.437)  # timestep in years

    # --- species 1: Ammophila ---
    #x1_spec1, x2_spec1, x3_spec1, x4_spec1, x5_spec1 = -0.6*dtYear, 0.2*dtYear, 0.6*dtYear, 2.0*dtYear, 2.2*dtYear
    #y1_spec1, y2_spec1, y3_spec1, y4_spec1, y5_spec1 = -0.5*dtYear, 0.0*dtYear, sp1_peak*dtYear, 0.0*dtYear, -0.5*dtYear
    
    # --- species 2: Elymus/Buckthorn ---
    x1_spec2, x2_spec2, x3_spec2, x4_spec2, x5_spec2 = -1.4*dtYear, -0.65*dtYear, 0.0*dtYear, 0.2*dtYear, 2.8*dtYear
    y1_spec2, y2_spec2, y3_spec2, y4_spec2, y5_spec2 = -0.5*dtYear, 0.0*dtYear, sp2_peak*dtYear, 0.0*dtYear, -0.5*dtYear

    # test with montlhy curves for ammophila
    # Ammophila (species 1)
    x1_spec1, y1_spec1 = -0.25, -0.30
    x2_spec1, y2_spec1 =  0, 0
    x3_spec1, y3_spec1 =  0.07,  0.12
    x4_spec1, y4_spec1 =  0.15,  0.02
    x5_spec1, y5_spec1 =  0.30, -0.10   
    
    # --- Sedimentation/erosion ---
    dh_df = 0.0
    if abs(srampnew-srampold) > 0.0:
        dh_df = srampnew-srampold
    elif abs(ysnew-ysold) > 0.0:
        dh_df = (ysnew-ysold)*betas

    dh_dt = snew - sold
    dh_db = -(ylnew-ylold)*betal

    # --- Slopes ---
    s12_spec1 = (y2_spec1-y1_spec1)/(x2_spec1-x1_spec1)
    s23_spec1 = (y3_spec1-y2_spec1)/(x3_spec1-x2_spec1)
    s34_spec1 = (y4_spec1-y3_spec1)/(x4_spec1-x3_spec1)
    s45_spec1 = (y5_spec1-y4_spec1)/(x5_spec1-x4_spec1)

    s12_spec2 = (y2_spec2-y1_spec2)/(x2_spec2-x1_spec2)
    s23_spec2 = (y3_spec2-y2_spec2)/(x3_spec2-x2_spec2)
    s34_spec2 = (y4_spec2-y3_spec2)/(x4_spec2-x3_spec2)
    s45_spec2 = (y5_spec2-y4_spec2)/(x5_spec2-x4_spec2)

    # --- Helper to compute growth increment with seasonal modulation ---
    def compute_increment(dh, x_vals, y_vals, s_vals, rhov_old):
        if dh < x_vals[0]:
            increment = 0.0
        elif dh >= x_vals[0] and dh < x_vals[1]:
            increment = (dh-x_vals[0])*s_vals[0] + y_vals[0]
        elif dh >= x_vals[1] and dh < x_vals[2]:
            increment = (dh-x_vals[1])*s_vals[1] + y_vals[1]
        elif dh >= x_vals[2] and dh < x_vals[3]:
            increment = (dh-x_vals[2])*s_vals[2] + y_vals[2]
        elif dh >= x_vals[3] and dh < x_vals[4]:
            increment = (dh-x_vals[3])*s_vals[3] + y_vals[3]
        else:
            increment = y_vals[4]

        # Apply seasonal factor only to positive increments
        if increment > 0 and month is not None:
            increment *= seasonal_factor(month)

        return max(rhov_old + increment, 0.0)  # ensure >= 0

    # --- Compute new vegetation ---
    rhov_dfnew_spec1 = compute_increment(dh_df, [x1_spec1, x2_spec1, x3_spec1, x4_spec1, x5_spec1],
                                         [y1_spec1, y2_spec1, y3_spec1, y4_spec1, y5_spec1],
                                         [s12_spec1, s23_spec1, s34_spec1, s45_spec1],
                                         rhov_df_spec1)
    rhov_dtnew_spec1 = compute_increment(dh_dt, [x1_spec1, x2_spec1, x3_spec1, x4_spec1, x5_spec1],
                                         [y1_spec1, y2_spec1, y3_spec1, y4_spec1, y5_spec1],
                                         [s12_spec1, s23_spec1, s34_spec1, s45_spec1],
                                         rhov_dt_spec1)
    rhov_dbnew_spec1 = compute_increment(dh_db, [x1_spec1, x2_spec1, x3_spec1, x4_spec1, x5_spec1],
                                         [y1_spec1, y2_spec1, y3_spec1, y4_spec1, y5_spec1],
                                         [s12_spec1, s23_spec1, s34_spec1, s45_spec1],
                                         rhov_db_spec1)

    # Species 2 only grows if species 1 is present
    if rhov_dfnew_spec1 > 0.5 or rhov_df_spec2 > 0.0:
        rhov_dfnew_spec2 = compute_increment(dh_df, [x1_spec2, x2_spec2, x3_spec2, x4_spec2, x5_spec2],
                                             [y1_spec2, y2_spec2, y3_spec2, y4_spec2, y5_spec2],
                                             [s12_spec2, s23_spec2, s34_spec2, s45_spec2],
                                             rhov_df_spec2)
    else:
        rhov_dfnew_spec2 = 0.0

    if rhov_dtnew_spec1 > 0.5 or rhov_dt_spec2 > 0.0:
        rhov_dtnew_spec2 = compute_increment(dh_dt, [x1_spec2, x2_spec2, x3_spec2, x4_spec2, x5_spec2],
                                             [y1_spec2, y2_spec2, y3_spec2, y4_spec2, y5_spec2],
                                             [s12_spec2, s23_spec2, s34_spec2, s45_spec2],
                                             rhov_dt_spec2)
    else:
        rhov_dtnew_spec2 = 0.0

    if rhov_dbnew_spec1 > 0.5 or rhov_db_spec2 > 0.0:
        rhov_dbnew_spec2 = compute_increment(dh_db, [x1_spec2, x2_spec2, x3_spec2, x4_spec2, x5_spec2],
                                             [y1_spec2, y2_spec2, y3_spec2, y4_spec2, y5_spec2],
                                             [s12_spec2, s23_spec2, s34_spec2, s45_spec2],
                                             rhov_db_spec2)
    else:
        rhov_dbnew_spec2 = 0.0

    # reduction of vegetation coverage due to frosting and reduction in winter
    
    if month is not None:
        # Define winter months (Dec, Jan, Feb)
        if month in [11, 12, 1, 2, 3]:
            winter_decay = 0.04  # reduction measuerd of 10-20% per winter, so about 
            winter_decay = 0.05  # reduction measuerd of 10-20% per winter, so about 

            rhov_dfnew_spec1 = max(0, rhov_dfnew_spec1 - winter_decay)
            rhov_dtnew_spec1 = max(0, rhov_dtnew_spec1 - winter_decay)
            rhov_dbnew_spec1 = max(0, rhov_dbnew_spec1 - winter_decay)
            rhov_dfnew_spec2 = max(0, rhov_dfnew_spec2 - winter_decay)
            rhov_dtnew_spec2 = max(0, rhov_dtnew_spec2 - winter_decay)
            rhov_dbnew_spec2 = max(0, rhov_dbnew_spec2 - winter_decay)
                       
            
    
    # Cap values between 0 and 1, and adjust species balance as before
    rhov_dfnew_spec1 = min(max(rhov_dfnew_spec1, 0.0), 1.0)
    rhov_dtnew_spec1 = min(max(rhov_dtnew_spec1, 0.0), 1.0)
    rhov_dbnew_spec1 = min(max(rhov_dbnew_spec1, 0.0), 1.0)
    rhov_dfnew_spec2 = min(max(rhov_dfnew_spec2, 0.0), 1.0)
    rhov_dtnew_spec2 = min(max(rhov_dtnew_spec2, 0.0), 1.0)
    rhov_dbnew_spec2 = min(max(rhov_dbnew_spec2, 0.0), 1.0)

    # Balance species
    if (rhov_dfnew_spec1 + rhov_dfnew_spec2 > 1.0):
        rhov_dfnew_spec1 = 1.0 - rhov_dfnew_spec2
    if (rhov_dtnew_spec1 + rhov_dtnew_spec2 > 1.0):
        rhov_dtnew_spec1 = 1.0 - rhov_dtnew_spec2
    if (rhov_dbnew_spec1 + rhov_dbnew_spec2 > 1.0):
        rhov_dbnew_spec1 = 1.0 - rhov_dbnew_spec2

    return rhov_dfnew_spec1, rhov_dtnew_spec1, rhov_dbnew_spec1, rhov_dfnew_spec2, rhov_dtnew_spec2, rhov_dbnew_spec2

def veg(ys, ysnew, yl, ylnew, sramp, srampnew, betas, betal, s, snew, dt,rhov_df,rhov_dt,rhov_db):
    """
    Vegetation growth function adapted from Duran and Moore (2014) 
    """

    HVmax = 1.0
    HVmin = 0.01 # treshold value
    tveg = 3*24*3600  # duran and moore [2013]
    tveg = 365*24*3600  # needs to be adjusted according to Bonte et al [2021]
    gamv = 0.05  # Duran and Moore[2013] default at 1.0 - plant sensitivity to erosion/accretion
    # start with ramp
    HV = max((rhov_df)**0.5*HVmax,HVmin) # vegetation height
    # change in bed level
    dh = 0.0
    if (abs(ysnew-ys)>0.0):
        dh=(ysnew-ys)*np.sin(betas)
    elif (abs(srampnew-sramp)>0.0):
        dh = 0.5*(srampnew-sramp)
        
#    rhov_dfnew=dt*(1-rhov_df)/tveg-gamv/HV*abs(dh)
    rhov_dfnew=rhov_df+dt*(1-rhov_df)/tveg-gamv*abs(dh)/HVmax
# Bonte
#    rhov_dfnew=rhov_df*1.005-0.005*rhov_df**2
    
    # update
    rhov_df = max(rhov_dfnew,0.0001)
    
    # follow up for the top of the dune
    HV = max((rhov_dt)**0.5*HVmax,HVmin) # vegetation height
    # change in bed level
    dh = (snew-s)
        
    rhov_dtnew=rhov_dt+dt*(1-rhov_dt)/tveg-gamv/HVmax*abs(dh)
    # update
    rhov_dt = max(rhov_dtnew,0.0001)
    
    # and the back of the dune
    HV = max((rhov_db)**0.5*HVmax,HVmin) # vegetation height
    # change in bed level
    dh=(ylnew-yl)*np.sin(betal)
    rhov_dbnew=rhov_db+dt*(1-rhov_db)/tveg-gamv/HVmax*abs(dh)
    # update
    rhov_db = max(rhov_dbnew,0.0001)
    
    return rhov_df, rhov_dt, rhov_db


### Morphodynamic functions
def sedbud2_1_nosmax(vramp, ys, ysprim, betas, s, yl, ylprim, betal, wtransp):
    # Dune evolution for positive sediment budget in case of triangular shape

    ylprim = ysprim #triangular shape, so landward and seaward crest are the same
    #snew = np.sqrt(((ys - yl) * s + 2 * wtransp) / (1 / betas + 1 / betal))
    snew = np.sqrt(s**2 + 2 * wtransp / (1 / betal + 1 / betas))
    ylprimnew = ylprim  # Dune increases in height but does not change in position
    ysprimnew = ysprim
    ylnew = ylprimnew - snew / betal
    ysnew = ysprimnew + snew / betas
    srampnew = snew
    vrampnew = 0.5 * (ysnew - ysprimnew)**2 * betas

    return snew, ylprimnew, ysprimnew, ylnew, ysnew, vrampnew, srampnew

def sedbud2_2_nosmax(vramp, ys, ysprim, betas, s, yl, ylprim, betal, wtransp):
    # dune evolution for positive sediment budget in case of trapezoidal shape
    vtrimax = 0.5 * (ysprim - ylprim) ** 2 / (1 / betal + 1 / betas)  # max volume to reach triangular shape
    vmax = vtrimax

    if wtransp <= vmax: #volume not enough to reach triangular shape
        snew = s + 2 * vtrimax / (ysprim - ylprim) - np.sqrt(2 * (vtrimax - wtransp) / (1 / betal + 1 / betas))
        srampnew = snew # do update here because you added the scarping part earlier
        #srampnew = sramp #do not update now because we are calculating in blocks
        ylnew = yl
        ylprimnew = ylprim + (snew - s) / betal
        ysprimnew = ysprim - (snew - s) / betas
        ysnew = ys
        vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas  #### do not recalculate vramp anymore
        #vrampnew = vramp  #keep the volume the same as we are calculating in blocks now

    else: #volume enough to reach triangular shape, transfer landward and seaward sediment
        wrest = wtransp - vmax

        #fill the triangular shape
        snew = s + 2 * vtrimax / (ysprim - ylprim)
        ylprimnew = ylprim + (snew - s) / betal
        ysprimnew = ysprim - (snew - s) / betas

        # expand the dune laterally using the remaining sediment
        #ysnew = ysprimnew + snew / betas
        #ylnew = ylprimnew - snew / betal
        #srampnew = snew
        #vrampnew = 0.5 * (ysnew - ysprimnew)**2 * betas

        # lee side
        ylnew = yl - wrest/snew 
        ylprimnew = ylnew + snew/betal

        # ramp side
        vrampnew = vramp + wrest/2 # get the new ramp volume

        # check the max volume for the ramp
        vrampmax = snew ** 2 / (2 * betas)

        if vrampnew < vrampmax: # check if the volume in the front of dune is below the maximum
            srampnew = (2 * betas * vrampnew)**0.5 # get the new ramp height --> should always be below s
            #ysnew = srampnew / betas + ysprim
            ysnew = ys #does not change because the ramp is only filling upwards
            ysprimnew = ysprim
            #vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas
        else: # if the ramp volume is larger than the maximum, move the ramp as a whole seaward, keeping the same angle
            extrasand = min(wrest/2, vrampnew - vrampmax) #get the excess of sand needed and redistribute it
            srampnew = snew # ramp is full, so it is equal to the crest height
            ysnew = (2 * extrasand / s) + ys # get the new seaward position of the ramp
            ysprimnew = ysnew - srampnew / betas # get the new crest position of the ramp 
            vrampnew = vrampmax # set the ramp volume to the maximum, because the extrasand becomes part of the core

    return snew, ylprimnew, ysprimnew, ylnew, ysnew, vrampnew, srampnew

def sedbud1_nosmax(vramp, ys, ysprim, betas, s, yl, ylprim, betal, wtransp):
    # dune evolution for stable sediment budget
    # print(f'in sebud1, wtransp', {wtransp})
    if ysprim > ylprim:  # TRAPEZOIDAL SHAPE, DISTRIBUTE SEDIMENT 50./50 ON CREST AND SEAWARD SLOPE
        # CHECK AMOUNT OF SEDIMENT REQUIRED TO ATTAIN TRIANGULAR SHAPE
        vtrimax = 0.5 * (ysprim - ylprim) ** 2 / (1 / betal + 1 / betas)

        if wtransp / 2 <= vtrimax:
            snew = s + (ysprim - ylprim) / (1 / betal + 1 / betas) - np.sqrt((ysprim - ylprim) ** 2 / (1 / betal + 1 / betas) ** 2 - wtransp / (1 / betas + 1 / betal))
            ylprimnew = ylprim + (snew - s) / betal
            ysprimnew = ysprim - (snew - s) / betas + 0.5 * wtransp / snew # just increase seaward
            ysnew = ysprimnew + snew / betas
            srampnew = snew
            ylnew = yl
            vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas

        elif wtransp / 2 > vtrimax:
            snew = s + 2 * vtrimax / (ysprim - ylprim) #get the maximum possible
            ylprimnew = ylprim + (snew - s) / betal #find the new position of the landward crest
            ysprimnew = ylprimnew # as triangular, seaward is the same as landward
            #ysnew = ysprimnew + snew / betas # added by filipe
            #ylnew = ylprimnew - snew / betal # added by filipe
            #srampnew = snew # added by filipe
            wrest = wtransp - vtrimax
            wtransp = wrest # amount of sediment left to distribute
            s = snew
            ysprimnew = ysprim  + 0.5 * wtransp / snew #added by filipe
            ylprimnew = ylprim  - 0.5 * wtransp / snew #added by filipe
            ysnew = ys + 0.5 * wtransp / snew #added by filipe so as to redistribute the leftovers on both sides
            ylnew = yl - 0.5 * wtransp / snew #added by filipe
            #ylprimnew = ylprim + (s / betal) # added by Filipe - recalculate after expanding the dune laterally
            #ysprimnew = ysprim - (s / betas)  # added by Filipe - recalculate after expanding the dune laterally
            ylprim = ylprimnew
            ysprim = ysprimnew
            vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas
            srampnew = (ysnew - ysprimnew) * betas  # added by Filipe
        # GOTO 1111

    elif ysprim <= ylprim + 0.001: # TRIANGULAR SHAPE
        snew = np.sqrt(((ys - yl) * s + 2 * wtransp) / (1 / betas + 1 / betal))
        ylprimnew = ylprim
        ysprimnew = ysprim
        ysnew = ysprimnew + snew / betas
        ylnew = ylprimnew - snew / betal
        srampnew = snew
        vrampnew = 0.5 * (ysnew - ysprimnew) ** 2 * betas

    return snew, ylprimnew, ysprimnew, ylnew, ysnew, vrampnew, srampnew

def sedbud0(ys, ysprim, betas, s, yl, ylprim, betal, wtransp):
    # dune evolution for negative sediment budget

    if ysprim > ylprim:
        vtrimax = 0.5 * (ysprim - ylprim) ** 2 * betal
        if wtransp / 2 <= vtrimax:
            ylprimnew = ysprim - np.sqrt((ysprim - ylprim) ** 2 - wtransp / betal)
            snew = s + (ylprimnew - ylprim) * betal
            ylprim = ylprimnew
            ylprimnew = ylprim - 0.5 * wtransp / snew
            ysprimnew = ysprim
            ysnew = ys
            ylnew = yl - 0.5 * wtransp / snew

        elif wtransp / 2 > vtrimax:
            snew = s + 2 * vtrimax / (ysprim - ylprim)
            ylprimnew = ylprim + (snew - s) / betal
            ysprimnew = ysprim
            wrest = wtransp - vtrimax
            wtransp = wrest
            s = snew
            snew = np.sqrt(betal * (s * (ysprim - yl) + 2 * wtransp))
            ysnew = ys
            ylnew = ylprimnew - snew / betal
    elif ysprim <= ylprim:
        snew = np.sqrt(betal * (s * (ysprim - yl) + 2 * wtransp))
        ysnew = ys
        ylprimnew = ylprim
        ysprimnew = ysprim
        ylnew = ylprimnew - snew / betal

    return snew, ylprimnew, ysprimnew, ylnew, ysnew

def duneeros(vramp, ysprim, betas, betal, s, erod, yl, ylprim, dfoot):
    # Calculate new ramp volume
    deltavramp = erod
    vrampnew = vramp + deltavramp

    if vrampnew > 0:
        srampnew = (2 * betas * vrampnew) ** 0.5
        ysnew = srampnew / betas + ysprim
        ysprimnew = ysprim
        ylprimnew = ylprim
        ylnew = yl
        snew = s

    elif vrampnew < 0:
        ysprimnew = ysprim + vrampnew / s

        if ysprimnew < ylprim:
            wtransprest = (ysprimnew - ylprim) * s
            ylprimnew = ((1 / betal) * ((ylprim - yl) * s + 2 * wtransprest)) ** 0.5 + yl
            snew = (ylprimnew - yl) * betal
            ysprimnew = ylprimnew
        else:
            ylprimnew = ylprim
            snew = s

        ysnew = ysprimnew
        ylnew = yl
        vrampnew = 0.0
        srampnew = 0.0

    else:  # vrampnew == 0
        ysnew = ysprim
        ysprimnew = ysprim
        ylprimnew = ylprim
        ylnew = yl
        srampnew = 0.0
        snew = s

    return vrampnew, ysnew, ysprimnew, snew, ylnew, ylprimnew, srampnew

### Transport functions
def wind_veg_reduction(u_in, L_block,
                   rhov1, Cd1=0.03,
                   rhov2=0.0, Cd2=0.01):
    """
    Wind speed reduction due to vegetation (can include two species).

    Parameters
    ----------
    u_in : float
        Incoming wind speed (m/s).
    L_block : float
        Length of the block (m).
    rhov1 : float
        Vegetation density (0–1) for species 1 (dune builder).
    Cd1 : float, optional
        Drag coefficient for species 1 (default=0.01).
    rhov2 : float, optional
        Vegetation density (0–1) for species 2 (stabilizer).
        Default = 0 (absent).
    Cd2 : float, optional
        Drag coefficient for species 2 (default=0.01).

    Returns
    -------
    u_out : float
        Outgoing wind speed after vegetation (m/s).
    """
    k = 0.5 * (Cd1 * rhov1 + Cd2 * rhov2)
    u_out = u_in * np.exp(-k * L_block)
    return u_out

def deposition_dune(dune_x, dune_z, substrate_x, substrate_z, qwnet, dt, tetaw, wspeed, rhos, rhoa, d50, kappa, z, vegcount, awind, rhov_df, rhov_dt, rhov_db, 
        rhov_df2, rhov_dt2, rhov_db2):
    """
    Calculate sediment deposition on dune in three blocks (ramp, crest, lee) based on wind shear velocity distribution.
    """

    dx_ramp = abs(dune_x['ramp_toe']-dune_x['ramp_crest']) # find the distance between the dunetoe and the dune crest
    dx_crest = abs(dune_x['ramp_crest']-dune_x['lee_crest']) #find the x distance between the front and back of the crest
    dx_lee = abs(dune_x['lee_foot']-dune_x['lee_crest']) # find the distance between the back of the crest and the dune lee foot
    
    # function to check where the dune is sitting
    def dike_surface_z(x, substrate_x, substrate_z):
        x_crest_start = substrate_x[2]
        x_crest_end   = substrate_x[3]
        x_toe         = substrate_x[4]
    
        z_crest = substrate_z[2]
        z_toe   = substrate_z[4]
    
        if x <= x_crest_start:
            # landward side (if relevant)
            return z_toe
    
        elif x_crest_start < x <= x_crest_end:
            # crest plateau
            return z_crest
    
        elif x_crest_end < x < x_toe:
            # seaward slope (linear face)
            return z_crest + (z_toe - z_crest) * (x - x_crest_end) / (x_toe - x_crest_end)
    
        else:
            # foreland
            return z_toe

    # check if everything should be deposited in the ramp or not
    z_dike = dike_surface_z(dune_x['seaward_crest'], substrate_x, substrate_z)
    
    if dx_ramp == 0 or dune_x['ramp_toe'] <= substrate_x[4] or (dune_x['seaward_crest'] <= substrate_x[3] and dune_z['seaward_crest'] <= substrate_z[3]) or abs(dune_z['seaward_crest'] - z_dike) <= 0.005 or dune_z['seaward_crest'] < z_dike - 0.005: # condition where no ramp exists anymore because of erosion, so not need to calculate everything as everything will end up on the ramp
        qs_ramp = qwnet * dt
        qs_crest = 0.0
        qs_lee = 0.0    
    else:
        
         # positions of each block center
        x_ramp  = (dune_x['ramp_toe'] + dune_x['ramp_crest']) / 2
        x_crest = (dune_x['ramp_crest'] + dune_x['lee_crest']) / 2
        x_lee   = (dune_x['lee_crest'] + dune_x['lee_foot']) / 2
        
         # compute slopes
        slope_ramp = dune_z['ramp_crest']/dx_ramp if dx_ramp != 0 else 0
        slope_crest = 0 # slope of the crest is zero as it is flat
        slope_lee = dune_z['seaward_crest']/dx_lee if dx_lee != 0 else 0
        
        
        # # block edges (note: keep these as signed positions along x)
        x_ramp_minus, x_ramp_plus   = dune_x['ramp_crest'], dune_x['ramp_toe']
        x_crest_minus, x_crest_plus = dune_x['lee_crest'], dune_x['ramp_crest']
        x_lee_minus, x_lee_plus     = dune_x['lee_foot'], dune_x['lee_crest']
         
        
         # ## exact per-block integrals (exclude self-block)
        I_ramp  = (slope_crest/np.pi) * np.log((x_ramp - x_crest_minus)/(x_ramp - x_crest_plus)) \
                 + (slope_lee  /np.pi) * np.log((x_ramp - x_lee_minus  )/(x_ramp - x_lee_plus   ))
         
        I_crest = (slope_ramp /np.pi) * np.log((x_crest - x_ramp_minus)/(x_crest - x_ramp_plus)) \
                 + (slope_lee  /np.pi) * np.log((x_crest - x_lee_minus )/(x_crest - x_lee_plus ))
         
        I_lee   = (slope_ramp /np.pi) * np.log((x_lee - x_ramp_minus)/(x_lee - x_ramp_plus)) \
                 + (slope_crest/np.pi) * np.log((x_lee - x_crest_minus)/(x_lee - x_crest_plus))
             
         # coefficients α (scaling factor controls amplitude of bed slope gradient influence) and β (relative wheigth of local slope term) for each block
        alpha_ramp, beta_ramp = 3, 0.8 #2, 0.8
        alpha_crest, beta_crest = 1.5, 0.2
        alpha_lee, beta_lee = 1.2, 0.8
        
         # calculate tau over tau0 including convolution term + local slope term + 1 (background shear)
        tau_ramp  = max(alpha_ramp * (I_ramp + beta_ramp * slope_ramp) + 1, 0.1)
        tau_crest = max(alpha_crest * (I_crest + beta_crest * slope_crest) + 1, 0.1)
        tau_lee   = max(alpha_lee * (I_lee + beta_lee * -slope_lee) + 1, 0.1)
                 
         # calculate the wind velocity distribution along the profile
        u10x_ramp = wspeed* np.cos(tetaw) * np.sqrt(tau_ramp)
        u10x_crest = wspeed * np.cos(tetaw) * np.sqrt(tau_crest)
        u10x_lee = wspeed * np.cos(tetaw) * np.sqrt(tau_lee)
        
        u10y = wspeed * np.sin(tetaw)
        
         # calculate critical velocities for sediment transport
        wcrit = awind * np.sqrt((rhos -rhoa)/rhoa * 9.81 * d50)
        tu10 = wcrit/ kappa * np.log(z / (d50/30))
        
         # allocate the bare velocities
        umag_ramp_bare = np.sqrt(u10x_ramp**2 + u10y**2)
        umag_crest_bare = np.sqrt(u10x_crest**2 + u10y**2)
        umag_lee_bare = np.sqrt(u10x_lee**2 + u10y**2)
        
         # Apply sequentially to x- and y-components
        if vegcount == 0:
            u10x_ramp  = wind_veg_reduction(u10x_ramp, dx_ramp, rhov_df)
            u10x_crest = wind_veg_reduction(u10x_crest, dx_crest, rhov_dt)
            u10x_lee   = wind_veg_reduction(u10x_lee, dx_lee, rhov_db)
        
            u10y_ramp  = wind_veg_reduction(u10y, dx_ramp, rhov_df)
            u10y_crest = wind_veg_reduction(u10y, dx_crest, rhov_dt)
            u10y_lee   = wind_veg_reduction(u10y, dx_lee, rhov_db)
        
        elif vegcount == 1:
            u10x_ramp  = wind_veg_reduction(u10x_ramp, dx_ramp, rhov_df, Cd1=0.12, rhov2=rhov_df2, Cd2=0.01)
            u10x_crest = wind_veg_reduction(u10x_crest, dx_crest, rhov_dt, Cd1=0.12, rhov2=rhov_dt2, Cd2=0.01)
            u10x_lee   = wind_veg_reduction(u10x_lee, dx_lee, rhov_db, Cd1=0.12, rhov2=rhov_db2, Cd2=0.01)
        
            u10y_ramp  = wind_veg_reduction(u10y, dx_ramp, rhov_df, Cd1=0.12, rhov2=rhov_df2, Cd2=0.01)
            u10y_crest = wind_veg_reduction(u10y, dx_crest, rhov_dt, Cd1=0.12, rhov2=rhov_dt2, Cd2=0.01)
            u10y_lee   = wind_veg_reduction(u10y, dx_lee, rhov_db, Cd1=0.12, rhov2=rhov_db2, Cd2=0.01)
        
        
         # Compute umag after reduction
        umag_ramp = np.sqrt(u10x_ramp**2 + u10y_ramp**2)
        umag_crest = np.sqrt(u10x_crest**2 + u10y_crest**2)
        umag_lee = np.sqrt(u10x_lee**2 + u10y_lee**2)
        umag_beach = wspeed * np.cos(tetaw) # calculation as in the subar routine
        

        # calculate sediment concentration potential
        cu_ramp = 1.5e-4 * (np.maximum(umag_ramp-tu10,0)**3) / np.maximum(umag_ramp, 1e-6)
        cu_crest = 1.5e-4 * (np.maximum(umag_crest-tu10,0)**3) / np.maximum(umag_crest, 1e-6)
        cu_lee = 1.5e-4 * (np.maximum(umag_lee-tu10,0)**3) / np.maximum(umag_lee, 1e-6)        
         
        cu_ramp_bare = 1.5e-4 * (np.maximum(umag_ramp_bare-tu10,0)**3) / np.maximum(umag_ramp_bare, 1e-6)
        cu_crest_bare = 1.5e-4 * (np.maximum(umag_crest_bare-tu10,0)**3) / np.maximum(umag_crest_bare, 1e-6)
        cu_lee_bare = 1.5e-4 * (np.maximum(umag_lee_bare-tu10,0)**3) / np.maximum(umag_lee_bare, 1e-6)
             
         # find the percentage of reduction due to vegetation
        F_c_ramp = 1 - min(cu_ramp / max(cu_ramp_bare, 1e-12), 1.0)
        F_c_crest = 1 -  min(cu_crest / max(cu_crest_bare, 1e-12), 1.0)
        F_c_lee = 1 - min(cu_lee / max(cu_lee_bare, 1e-12), 1.0)
         
        #find gradients
        Fg_ramp  = max(0, tu10 - umag_ramp) / max(tu10, 1e-12)
        Fg_crest = max(0, umag_ramp - umag_crest) / max(umag_ramp, 1e-12)
        Fg_lee   = max(0, umag_crest - umag_lee) / max(umag_crest, 1e-12)

        #Fg_ramp  = max(0, cu_beach - cu_ramp) / max(cu_beach, 1e-12)
        #Fg_crest = max(0, cu_ramp - cu_crest) / max(cu_ramp, 1e-12)
        #Fg_lee   = max(0, cu_crest - cu_lee) / max(cu_crest, 1e-12)
        
        Fg_sum = Fg_ramp + Fg_crest + Fg_lee
        
        #combined prob
        Fstar_ramp  = 1 - (1 - F_c_ramp)  * (1 - (Fg_ramp/max(Fg_sum, 1e-12)))
        Fstar_crest = 1 - (1 - F_c_crest) * (1 - (Fg_crest/max(Fg_sum, 1e-12)))
        Fstar_lee   = 1 - (1 - F_c_lee)   * (1 - (Fg_lee/max(Fg_sum, 1e-12)))

        qs_total = qwnet * dt

        qs_ramp  = qs_total * Fstar_ramp
        rem      = qs_total - qs_ramp
        
        qs_crest = rem * Fstar_crest
        qs_lee   = qs_total - qs_ramp - qs_crest

         # apply gradient based, so to ensure deposition accross when bare
        # G_ramp  = max(0, tu10 - umag_ramp)
        # G_crest = max(0, umag_ramp  - umag_crest)
        # G_lee   = max(0, umag_crest - umag_lee)
         
        # G_sum = G_ramp + G_crest + G_lee
         
        # qs_total = qwnet * dt
         
        #  # ramp 
        # if F_c_ramp > 1e-6:  # vegetation present → use it
        #     qs_ramp = qs_total * F_c_ramp
        # else:                # no veg → fallback to gradient
        #     qs_ramp = qs_total * (G_ramp / max(G_sum, 1e-12))
         
        # # --- crest ---
        # rem_after_ramp = qs_total - qs_ramp
        # if F_c_crest > 1e-6:
        #     qs_crest = rem_after_ramp * F_c_crest
        # else:
        #     qs_crest = rem_after_ramp * (G_crest / max(G_crest + G_lee, 1e-12))
         
        # # --- lee (whatever remains) ---
        # qs_lee = qs_total - (qs_ramp + qs_crest)
        
         
        if dune_z['ramp_crest'] < dune_z['seaward_crest']: #case with scarp, enforce all depo in front             
            qs_ramp = qwnet * dt
            qs_crest = 0.0
            qs_lee = 0.0

    return qs_ramp, qs_crest, qs_lee

def morpho_update(dune_x, dune_z, substrate_x, substrate_z, qs_ramp, qs_crest, qs_lee, qs, ql, r, swl, vnourd, dt,
                   betas, betas_hs, betal, betal_hs, vramp, deltavall):
    """
    Update dune and substrate morphology based on sediment deposition and erosion.
    """

    erod = (qs + ql) * dt
    erod_ramp = qs * dt
    erod_overwash = ql * dt
    
    # get dune and dike parameters that to enter the calculations
    s_hs = substrate_z[3]-substrate_z[4]  # dike height above substrate seafloor
    ylprim_hs = substrate_x[2]  # landward dike crest position
    ysprim_hs = substrate_x[3]  # seaward dike crest position
    ys_hs = substrate_x[4]  # dike seaward toe position
    dfoot = dune_z['ramp_toe']  # height of dune foot above substrate seafloor
    ys = dune_x['ramp_toe'] # dune seaward toe position
    yl = dune_x['lee_foot']  # dune landward toe position
    s = dune_z['seaward_crest']-dfoot  # dune height above
    sramp = dune_z['ramp_crest']-dfoot  # dune ramp height above
    ylprim = dune_x['lee_crest']  # dune landward crest position
    ysprim = dune_x['seaward_crest']  # dune seaward crest position


    # case 1: dune erosion
    # case 1: dune erosion
    if erod > 0:
        
        srampnew = dune_z['ramp_crest']-dfoot #allocate it if it does not change for vramp calc

        #check if there is dune to be eroded:
        if r + swl <= (s_hs+dfoot): # if water level between dike foot and dike crest
            z_seg = substrate_z[3:5]
            x_seg = substrate_x[3:5]
            # sort by z
            order = np.argsort(z_seg)
            z_sorted = np.array(z_seg)[order]
            x_sorted = np.array(x_seg)[order]
            runup_reach_dike_x = np.interp(r + swl, z_sorted, x_sorted)

            z_seg = [dune_z['ramp_toe'], dune_z['ramp_crest']]
            x_seg = [dune_x['ramp_toe'], dune_x['ramp_crest']]
            # sort by z
            order = np.argsort(z_seg)
            z_sorted = np.array(z_seg)[order]
            x_sorted = np.array(x_seg)[order]
            runup_reach_dune_x = np.interp(r + swl, z_sorted, x_sorted)
            
            if runup_reach_dike_x >= runup_reach_dune_x: #dike exposed
                erod_ramp = 0
            
            if erod_ramp > 0: #if there is potential erosion on the dune ramp
                # calculate maximum erosion volume on the ramp based on the dike limit
                volume_max_ero = partial_volume(dune_x, dune_z, substrate_x, substrate_z, runup_reach_dike_x, dune_x['ramp_toe'])
                if erod_ramp >= volume_max_ero:
                    erod_ramp = volume_max_ero #pad it to the max erosion volume on the ramp
                
                def f_ramp_erosion(x):
                    volume_potential = partial_volume(dune_x, dune_z, substrate_x, substrate_z, x, dune_x['ramp_toe'])
                    return volume_potential - erod_ramp
                
                #find new locations that satisfy the erosion volume
                ysprim_new = brentq(f_ramp_erosion, runup_reach_dike_x, dune_x['ramp_toe']) #brettq to find the new dune crest position after erosion
                #if ysprim_new <= runup_reach_dike_x: 
                srampnew = max(np.interp(ysprim_new, list(dune_x.values()), list(dune_z.values())), np.interp(ysprim_new, substrate_x, substrate_z))
                #else:
                #    srampnew = dune_z['ramp_toe'] #no change as the dike is exposed

        else: #dike submerged thus overtopping as well
            
            #if ysprim_hs >= dune_x['ramp_crest'] and dune_z['seaward_crest'] - dfoot <= s_hs:
            #    erod_ramp = 0

            if dune_x['ramp_crest'] <= ylprim_hs and dune_z['seaward_crest'] - dfoot > s_hs:
                erod_overwash = erod_ramp + erod_overwash #assume all go behind because there is no return flow anymore
                erod_ramp = 0
            if erod_ramp > 0: #if there is potential erosion on the dune ramp
            # calculate maximum erosion volume on the ramp based on the dike limit
                volume_max_ero = partial_volume(dune_x, dune_z, substrate_x, substrate_z, dune_x['lee_foot'], dune_x['ramp_toe']) #whole dune is game
                
                if erod_ramp >= volume_max_ero:
                    erod_ramp = volume_max_ero #pad it to the max erosion volume on the ramp

                def f_ramp_erosion(x):
                    volume_potential = partial_volume(dune_x, dune_z, substrate_x, substrate_z, x, dune_x['ramp_toe'])
                    return volume_potential - erod_ramp
                
                #find new locations that satisfy the erosion volume
                ysprim_new = brentq(f_ramp_erosion, dune_x['lee_foot'], dune_x['ramp_toe']) #brettq to find the new dune crest position after erosion
                #if ysprim_new <= ys_hs:
                srampnew = max(np.interp(ysprim_new, list(dune_x.values()), list(dune_z.values())), np.interp(ysprim_new, substrate_x, substrate_z))
                #else:
                #    srampnew = dune_z['ramp_toe'] #no change as the dike is exposed

        deltavramp = erod_ramp
        deltavall=deltavall+deltavramp

        if erod_ramp > 0:
        # update dune profile after erosion
            if ysprim_new > dune_x['ramp_crest'] and ysprim_new > ys_hs: #still old ramp to erode further, so it becomes the toe. The case for outside the dike is to enforce proper enter only without interaction. 
                dune_x['ramp_toe'] = ysprim_new
                dune_x['ramp_toe_temp'] = ysprim_new
                dune_z['ramp_toe_temp'] = srampnew
            elif ysprim_new <= dune_x['ramp_crest'] and ysprim_new > dune_x['lee_crest'] and ysprim_new > ys_hs: #erosion beyond the ramp but not beyond the dike
                dune_x['ramp_toe'] = ysprim_new
                dune_x['ramp_toe_temp'] = ysprim_new
                dune_x['ramp_crest'] = ysprim_new
                dune_z['ramp_crest'] = srampnew
                dune_x['ramp_dike_temp'] = ysprim_new
                dune_z['ramp_dike_temp'] = srampnew
                dune_x['seaward_crest'] = ysprim_new
                dune_x['crest_temp'] = ysprim_new
            elif ysprim_new <= dune_x['ramp_crest'] and ysprim_new > dune_x['lee_crest'] and ysprim_new <= ys_hs: #erosion beyond the ramp and within the dike
                dune_x['ramp_toe'] = ys_hs
                if ysprim_new > ysprim_hs:
                    dune_x['ramp_toe_temp'] = ysprim_new
                    dune_x['ramp_crest'] = ysprim_new
                    dune_z['ramp_crest'] = srampnew
                else:
                    dune_x['ramp_toe_temp'] = ysprim_hs
                    dune_x['ramp_crest'] = ysprim_new
                    dune_z['ramp_crest'] = np.interp(ysprim_new, substrate_x, substrate_z)
                dune_z['ramp_toe_temp'] = np.interp(ysprim_new, substrate_x, substrate_z)
                dune_x['ramp_dike_temp'] = ysprim_new
                dune_z['ramp_dike_temp'] = srampnew
                dune_x['seaward_crest'] = ysprim_new
                dune_x['crest_temp'] = ysprim_new
                if ysprim_new < ylprim_hs:
                    dune_x['ramp_toe'] = ys_hs
                    dune_z['ramp_toe'] = np.interp(ys_hs, substrate_x, substrate_z)
                    dune_x['ramp_toe_temp'] = ysprim_hs
                    dune_z['ramp_toe_temp'] = np.interp(ysprim_hs, substrate_x, substrate_z)
                    dune_x['ramp_crest'] = ylprim_hs
                    dune_z['ramp_crest'] = np.interp(ylprim_hs, substrate_x, substrate_z)
                    dune_x['ramp_dike_temp'] = ysprim_new
                    dune_z['ramp_dike_temp'] = np.interp(ysprim_new, substrate_x, substrate_z)

            elif ysprim_new <= dune_x['lee_crest'] and ysprim_new <= ys_hs: #erosion beyond the dune
                dune_x['ramp_toe'] = ysprim_new
                dune_z['ramp_toe'] = srampnew
                dune_x['ramp_toe_temp'] = ysprim_new
                dune_z['ramp_toe_temp'] = srampnew
                dune_x['ramp_crest'] = ysprim_new
                dune_z['ramp_crest'] = srampnew
                dune_x['ramp_crest_temp'] = ysprim_new
                dune_z['ramp_crest_temp'] = srampnew
                dune_x['ramp_dike_temp'] = ysprim_new
                dune_z['ramp_dike_temp'] = srampnew
                dune_x['lee_crest'] = ysprim_new
                dune_z['lee_crest'] = srampnew
                dune_x['crest_temp'] = ysprim_new
                dune_z['crest_temp'] = srampnew
                dune_x['seaward_crest'] = ysprim_new
                dune_z['seaward_crest'] = srampnew

        deltavramp_ero = deltavramp
        vramp = partial_volume(dune_x, dune_z, substrate_x, substrate_z, dune_x['ramp_crest'], dune_x['ramp_toe'])
        
        vdune_test_before_overwash = partial_volume(dune_x, dune_z, substrate_x, substrate_z, dune_x['lee_foot'], dune_x['ramp_toe'])
        
        if erod_overwash > 0: # if some of this erosion is due to overwash, move this volume to the lee portion --> No change in total dune volume, only reorganization within the dune
            # check available crest volume to be transported to the lee side
            def f_crest_erosion(s_potential):
                # Volume of the triangular cap removed when lowering the crest by s_potential
                delta = s_potential
                crest_width = dune_x['seaward_crest'] - dune_x['lee_crest']
                V_cap = crest_width * delta + delta**2 / 2 * (1 / betas + 1 / betal)
                
                return V_cap - erod_overwash
            
            q_to_lee = 0.
                            
            s_new_ero = dune_z['seaward_crest'] - brentq(f_crest_erosion, 0, dune_z['seaward_crest'] - dfoot)
        
            #check if there is a dike and, if so, limit the erosion until the top of the dike
            if s_hs <= dune_z['seaward_crest'] - dfoot and ((ysprim_hs < dune_x['seaward_crest'] and ysprim_hs > dune_x['lee_crest']) or (ylprim_hs < dune_x['seaward_crest'] and ylprim_hs > dune_x['lee_crest'])) and s_new_ero - dfoot <= s_hs:
                if dune_z['seaward_crest'] - dfoot == s_hs:
                    q_to_lee = 0
                else:
                    #q_to_lee = f_crest_erosion(s_hs+dfoot-s_new_ero) + erod_overwash
                    q_to_lee = f_crest_erosion(dune_z['seaward_crest'] - (s_hs + dfoot)) + erod_overwash
                erod_overwash = q_to_lee
                s_new_ero = s_hs + dfoot #on these condition you cannot erode the dike
                ysprim_new_c = dune_x['seaward_crest'] - (dune_z['seaward_crest'] - s_new_ero) / betas
                ylprim_new_c = dune_x['lee_crest'] - (dune_z['seaward_crest'] - s_new_ero) / betal


                if dune_x['seaward_crest'] == dune_x['ramp_dike_temp'] and dune_x['seaward_crest'] >= ylprim_hs:
                    dune_x['ramp_dike_temp'] = ysprim_new_c
                    dune_z['ramp_dike_temp'] = s_new_ero
                    
                dune_x['seaward_crest'] = ysprim_new_c
                dune_x['lee_crest'] = ylprim_new_c
                dune_z['seaward_crest'] = s_new_ero
                dune_z['lee_crest'] = s_new_ero
                dune_z['crest_temp'] = s_new_ero
                dune_x['crest_temp'] = ysprim_new_c
                
                
                
            else:
                q_to_lee = erod_overwash
                ysprim_new_c = dune_x['seaward_crest'] - (dune_z['seaward_crest']- s_new_ero) / betas
                ylprim_new_c = dune_x['lee_crest'] - (dune_z['seaward_crest'] - s_new_ero) / betal
                
                if dune_x['seaward_crest'] == dune_x['ramp_dike_temp'] and dune_x['seaward_crest'] >= ylprim_hs:
                    dune_x['ramp_dike_temp'] = ysprim_new_c
                    dune_z['ramp_dike_temp'] = s_new_ero

                if dune_z['ramp_crest'] == dune_z['seaward_crest']:
                    dune_z['ramp_crest'] = s_new_ero

                if dune_x['seaward_crest'] < ylprim_hs: # if the dune crest is now behind the dike crest, adjust to the dike crest
                    #dune_z['ramp_dike_temp'] = dune_z['seaward_crest']
                    #dune_x['ramp_dike_temp'] = np.interp(dune_z['seaward_crest'], substrate_z[1:3], substrate_x[1:3])
                    dune_z['ramp_dike_temp'] = s_new_ero
                    dune_x['ramp_dike_temp'] = np.interp(s_new_ero, substrate_z[1:3], substrate_x[1:3])
                    
                dune_x['seaward_crest'] = ysprim_new_c
                dune_x['lee_crest'] = ylprim_new_c
                dune_z['seaward_crest'] = s_new_ero
                dune_z['lee_crest'] = s_new_ero
                dune_z['crest_temp'] = s_new_ero
                dune_x['crest_temp'] = ysprim_new_c

                    
            if q_to_lee > 0:
                deltavlee = q_to_lee
                
                ################## update lee ##################

                vlee_c = lee_volume(dune_x['lee_crest'], dune_x, dune_z, betal, substrate_x, substrate_z)
                vlee_target = vlee_c + deltavlee
                
                    
                def f_lee_overwash(ylprim_potential):
                    v_over = lee_volume(ylprim_potential, dune_x, dune_z, betal, substrate_x, substrate_z)
                    
                    return v_over - vlee_target
                
                # solve for ylprim_potential that gives the target volume
                ylprim_new_l = brentq(f_lee_overwash, dune_x['lee_crest'] - 100, dune_x['lee_crest']) # find the new lee height using the root finder
                yl_flat = ylprim_new_l - (dune_z['seaward_crest']-dune_z['lee_foot']) / betal # get the new lee foot position assuming no hard structure
                
                # use shapely to find intersection between new lee line and potential hard structure
                substrate_line = LineString(zip(substrate_x, substrate_z))
                # lee line: from lee crest (ylprim, s_new) to theoretical flat toe (yl_flat, dfoot)
                lee_line = LineString([(ylprim_new_l, dune_z['seaward_crest']), (yl_flat, dune_z['lee_foot'])])
                
                # find intersection(s)
                intersection = substrate_line.intersection(lee_line)
                
                if not intersection.is_empty:
                    if intersection.geom_type == "Point":
                        yl_new_l, z_int_l = intersection.x, intersection.y
                    else:  # multiple intersections, take the shallowest one (max x)
                        points = [pt for pt in intersection.geoms if pt.geom_type == "Point"]
                        yl_new_l, z_int_l = max(points, key=lambda p: p.x).x, max(points, key=lambda p: p.x).y
                
                    
                else:
                    # no intersection with hard substrate → use flat geometry
                    yl_new_l = yl_flat
                    z_int_l = dune_z['lee_foot'] 
                
                # Update dune profile
                # dune_x[1:3] = [yl_new_l, ylprim_new_l]
                # dune_z[1:3] = [z_int_l, s_new_c + dfoot]
                dune_x['lee_foot'] = yl_new_l
                dune_z['lee_foot'] = z_int_l
                dune_x['lee_crest'] = ylprim_new_l
                dune_z['lee_crest'] = s_new_ero

        vdune_test_after_overwash = partial_volume(dune_x, dune_z, substrate_x, substrate_z, dune_x['lee_foot'], dune_x['ramp_toe'])
        
        print('Delta Volume Overwash = ', vdune_test_after_overwash - vdune_test_before_overwash)
    # case 2: dune accretion
    
    else:
        
        deltavramp_ero = 0.0
        deltavramp = 0.0
            
        # Initialize the idx for volume calculations
        ramp_idx = ['ramp_crest', 'ramp_toe']
        crest_idx = ['lee_crest', 'crest_temp', 'seaward_crest']
        lee_idx = ['lee_foot', 'lee_crest']

        # chech dune nourishment
        if vnourd > 0:
            qs_ramp = qs_ramp + vnourd
            
        # change in the ramp
        if qs_ramp > 0:
            deltavramp = qs_ramp # how much sand is deposited on the ramp
            
            # current, max, new
            vramp      = ramp_volume(dune_z['ramp_crest']-dfoot, dune_z['ramp_crest']-dfoot, dune_x['ramp_crest'], dune_x['ramp_toe'], betas, dfoot, substrate_x, substrate_z)
            x_foot_max = dune_x['seaward_crest'] + ((dune_z['seaward_crest']-dfoot)/betas)
            vrampmax   = ramp_volume(dune_z['seaward_crest']-dfoot, dune_z['seaward_crest']-dfoot, dune_x['seaward_crest'], x_foot_max, betas, dfoot, substrate_x, substrate_z)
            vrampnew   = vramp + deltavramp
            
            ys_temp = dune_x['seaward_crest'] + ((dune_z['seaward_crest']-dfoot) /betas)
                    
            # check if ramp is in front of the dike or not:
            if (dune_x['ramp_crest'] > dune_x['seaward_crest']) or (dune_x['seaward_crest'] == dune_x['ramp_crest'] and dune_x['seaward_crest'] < ys_hs):  # if there is a ramp building in front of a exposed dike (hence to connections) or a flat surface connecting somewhere in the dike with ys = ysprim
                
                if dune_z['seaward_crest'] <= s_hs + dfoot: # there is no active dune on top of the dike, and everything is exposed
                    vrampmax_hs = ((((s_hs)**2)/(2*betas)))-((((s_hs)**2)/(2*betas_hs))) # volume to fill the front of the dike
                    x_rampfoot_max = substrate_x[3] + (s_hs)/betas
                    
                    v_dike = (substrate_x[4]-substrate_x[3])*(substrate_z[3]-substrate_z[4])/2
                    v_bottom_full = (x_rampfoot_max - substrate_x[3])*(substrate_z[3]-substrate_z[4])/2
                    v_bottom = ((dune_x['ramp_toe'] - dune_x['ramp_crest'])*(dune_z['ramp_crest'] - dune_z['ramp_toe'])/2) - ((substrate_x[4]- dune_x['ramp_crest'])*(dune_z['ramp_crest'] - dune_z['ramp_toe'])/2)
                    
                    if v_bottom < 0:
                        v_bottom = 0
                        
                    v_to_fill = v_bottom_full - v_dike - v_bottom #this is the volume to fill the front of the dike
                    
                    if deltavramp < v_to_fill: # if the volume is not sufficient to fill the whole front of the dike:
                        sramp_new_r = ((2 * (deltavramp+v_bottom)) / ((1 / betas) - (1 / betas_hs))) ** 0.5
                        ysprim_new_r = ys_hs - (sramp_new_r/betas_hs)
                        ys_new_r = ysprim_new_r+(sramp_new_r/betas)
                        
                        # check where to update:
                        if ((dune_z['seaward_crest'] - dune_z['ramp_toe']) > s_hs and dune_x['ramp_toe'] <= substrate_x[4]) or dune_x['seaward_crest'] < dune_x['ramp_crest']: #check if no ramp or detacehd
                           
                           #check if there is ramp on top of the dune or not
                           v_temp = partial_volume(dune_x, dune_z, substrate_x, substrate_z, dune_x['seaward_crest'], substrate_x[4])
                           
                           if v_temp > 0: #means there is a ramp on top, not scarp
                               dune_x['ramp_dike_temp'] = dune_x['seaward_crest'] + ((dune_z['seaward_crest'] - substrate_z[3])/betas)
                               dune_z['ramp_dike_temp'] = s_hs + dfoot
                           else: #assume a scarp
                               dune_x['ramp_dike_temp'] = dune_x['seaward_crest']
                               dune_z['ramp_dike_temp'] = s_hs + dfoot
       
                           dune_x['ramp_toe'] = ys_new_r
                           dune_x['ramp_toe_temp'] = ys_new_r
                           dune_z['ramp_crest'] = sramp_new_r+dfoot
                           dune_x['ramp_crest'] = ysprim_new_r

                        else:
                           
                            if dune_z['seaward_crest'] < substrate_z[3]: #check again if there is on top or not
                                dune_z['seaward_crest'] = sramp_new_r+dfoot
                                dune_z['ramp_dike_temp'] = sramp_new_r+dfoot
                                
                            dune_x['ramp_toe'] = ys_new_r
                            dune_x['ramp_toe_temp'] = ys_new_r
                            dune_z['ramp_crest'] = sramp_new_r+dfoot
                            dune_x['ramp_crest'] = ysprim_new_r
                            dune_x['seaward_crest'] = ysprim_new_r
                            dune_x['ramp_dike_temp'] = ysprim_new_r
                           
                        
                    elif v_to_fill < 0: #case where the front is already filled but enters this part because of the remaining pivot crest locations
                        #since there is no crest, just move the ramp forward
                        vextra = deltavramp # extra volume that needs to be accomodated by moving the dune toe and ysprim
                        
                        ys_new_r = (vextra / dune_z['ramp_crest']) + ys # get the new seaward position of the ramp
                        ysprim_new_r = ysprim + (vextra / dune_z['ramp_crest'])  # get the new crest position of the ramp 
                        
                        dune_x['ramp_toe'] = ys_new_r
                        dune_x['ramp_crest'] = ysprim_new_r
                        dune_x['seaward_crest'] = ysprim_new_r
                        dune_x['ramp_dike_temp'] = ysprim_new_r
                        dune_z['seaward_crest'] = dune_z['ramp_crest']
                        dune_z['ramp_dike_temp'] = dune_z['ramp_crest']
                        
                    else:
                        vextra = deltavramp-v_to_fill #amount to be used to redistribute it in front of the dike
                        sramp_new_r = s_hs #fill completely up to the dike height                
                        ysprim_new_r = ysprim + (vextra / sramp_new_r)  # get the new crest position of the ramp 
                        ys_new_r = (sramp_new_r/betas) + ysprim_new_r
                        
                        dune_x['ramp_toe'] = ys_new_r
                        dune_z['ramp_crest'] = sramp_new_r+dfoot
                        dune_x['ramp_crest'] = ysprim_new_r
                        dune_x['seaward_crest'] = ysprim_new_r
                        dune_z['seaward_crest'] = sramp_new_r+dfoot
                        dune_x['ramp_dike_temp']= ysprim_new_r #ysprim
                        dune_z['ramp_dike_temp'] = sramp_new_r+dfoot #s+dfoot
                        
                elif ys_temp < ys_hs and dune_z['seaward_crest'] > s_hs + 3: #there is dune on top but the height is not enough to fill the front
                    
                    # check if the volume needed to reconnect the ramp side

                    vrampmax_hs = ((((s_hs)**2)/(2*betas)))-((((s_hs)**2)/(2*betas_hs))) # volume to fill the front of the dike
                    x_rampfoot_max = substrate_x[3] + (s_hs)/betas
                    
                    v_dike = (substrate_x[4]-substrate_x[3])*(substrate_z[3]-substrate_z[4])/2
                    v_bottom_full = (x_rampfoot_max - substrate_x[3])*(substrate_z[3]-substrate_z[4])/2
                    v_bottom = ((dune_x['ramp_toe'] - dune_x['ramp_crest'])*(dune_z['ramp_crest'] - dune_z['ramp_toe'])/2) - ((substrate_x[4]- dune_x['ramp_crest'])*(dune_z['ramp_crest'] - dune_z['ramp_toe'])/2)
                    
                    if v_bottom < 0:
                        v_bottom = 0
                        
                    v_to_fill = v_bottom_full - v_dike - v_bottom #this is the volume to fill the front of the dike
                    
                    if deltavramp < v_to_fill: # if the volume is not sufficient to fill the whole front of the dike:
                       sramp_new_r = ((2 * (deltavramp+v_bottom)) / ((1 / betas) - (1 / betas_hs))) ** 0.5
                       ysprim_new_r = ys_hs - (sramp_new_r/betas_hs)
                       ys_new_r = ysprim_new_r+(sramp_new_r/betas)
                       
                       # check where to update:
                       if ((dune_z['seaward_crest'] - dune_z['ramp_toe']) > s_hs and dune_x['ramp_toe'] <= substrate_x[4]) or dune_x['seaward_crest'] < dune_x['ramp_crest']: #check if no ramp or detacehd
                           
                           #check if there is ramp on top of the dune or not
                           v_temp = partial_volume(dune_x, dune_z, substrate_x, substrate_z, dune_x['seaward_crest'], substrate_x[4])
                           
                           if v_temp > 0: #means there is a ramp on top, not scarp
                               dune_x['ramp_dike_temp'] = dune_x['seaward_crest'] + ((dune_z['seaward_crest'] - substrate_z[3])/betas)
                               dune_z['ramp_dike_temp'] = s_hs + dfoot
                           else: #assume a scarp
                               dune_x['ramp_dike_temp'] = dune_x['seaward_crest']
                               dune_z['ramp_dike_temp'] = s_hs + dfoot
       
                           dune_x['ramp_toe'] = ys_new_r
                           dune_x['ramp_toe_temp'] = ys_new_r
                           dune_z['ramp_crest'] = sramp_new_r+dfoot
                           dune_x['ramp_crest'] = ysprim_new_r

                       else:
                           dune_x['ramp_toe'] = ys_new_r
                           dune_x['ramp_toe_temp'] = ys_new_r
                           dune_z['ramp_crest'] = sramp_new_r+dfoot
                           dune_x['ramp_crest'] = ysprim_new_r
                           dune_x['seaward_crest'] = ysprim_new_r
                           dune_x['ramp_dike_temp'] = ysprim_new_r
                           
                    else:
                        #vol_to_fill_top = ((ysprim_hs - dune_x['seaward_crest'])**2)*betas*0.5 #volume to fill the top of the dike and stop the intersection
                        vol_to_fill_top = (ysprim_hs - (2*(dune_z['seaward_crest']-dune_z['ramp_crest'])/betas))*((dune_z['seaward_crest']-dune_z['ramp_crest'])/2)
                        vol_to_fill_face = ((((s_hs)**2)/(2*betas)))-((ys_hs-ysprim_hs)*(s_hs)*0.5) #volume to fill the face of the dike
                        vol_extra = vrampnew - vol_to_fill_face #extra volume that needs to be accomodated above the dike
                        
                            
                        if vol_extra < vol_to_fill_top: #thus not enough to fill everything and stop connection with dike
        
                            #check if the front is already filled, so just move all to the top (case 9b)
                            if v_to_fill < 0:
                                vol_extra = deltavramp
                                
                            def f_ramp_dike_top(sramp_potential):
                                #toebound = ysprim_hs
                                toebound = dune_x['ramp_crest'] # need to be the crest because once it starts growing, the volume limits
                                crestbound = dune_x['seaward_crest']
                                sramp_prim_x = toebound - (sramp_potential/betas)
                                dune_x_crest = np.array([crestbound, sramp_prim_x, toebound])
                                dune_z_crest = np.array([sramp_potential, sramp_potential, 0])
        
                                # volume between crest line and substrate
                                V = np.trapz(dune_z_crest, dune_x_crest)
                        
                                return V - vol_extra
                            
                            # solve for sramp_potential that gives the target volume
                            sramp_new_r = brentq(f_ramp_dike_top, 0, dune_z['seaward_crest']-dune_z['ramp_crest']) + dune_z['ramp_crest'] # find the new ramp height using the root finder padding only the top of the dike
                            dune_x['ramp_dike_temp'] = ysprim + (((dune_z['seaward_crest'])-sramp_new_r)/betas) #pad the previous so that it only updates when fully filled
                            dune_x['ramp_crest'] = dune_x['ramp_crest'] - ((sramp_new_r-dune_z['ramp_crest'])/betas) #foot based on the dike because it only moves once everything is filled
                            dune_z['ramp_crest'] = sramp_new_r
                            #dune_x['ramp_toe'] = dune_x['ramp_crest']+(sramp_new_r/betas)
                            #dune_z['ramp_toe'] = dfoot
                            dune_z['ramp_dike_temp'] = sramp_new_r #get the same height as the ramp until filled
        
                            
                        else: #enough to fill everything and stop connection with dike
                        
                            #find volume extra to be redistributed in the crest:
                            vol_extra = vol_extra-vol_to_fill_top
                            
                            # relink everything
                            dune_z['ramp_dike_temp'] = dune_z['seaward_crest']
                            dune_z['ramp_crest'] = dune_z['seaward_crest']
                            dune_x['ramp_crest'] = dune_x['ramp_toe'] - ((dune_z['ramp_crest']-dune_z['ramp_toe'])/betas)
                            dune_x['ramp_dike_temp'] =  dune_x['ramp_crest']
                            dune_x['seaward_crest'] =  dune_x['ramp_crest']
                            
                            #use the remaining volume to fill the crest
                            qs_crest = qs_crest+vol_extra
                    
                elif vrampnew < vrampmax: # there is some dune on the top or partial front
                    def f_ramp_dike(sramp_potential):
                        return ramp_volume(sramp_potential, dune_z['ramp_crest']-dfoot, dune_x['ramp_crest'], dune_x['ramp_toe'], betas, dfoot, substrate_x, substrate_z) - vrampnew
                    
                    # solve for sramp_potential that gives the target volume
                    sramp_new_r = brentq(f_ramp_dike, dune_z['ramp_crest']-dfoot, dune_z['seaward_crest']-dfoot) # find the new ramp height using the root finder
        
                    if sramp_new_r < s_hs: #check if the sramp is within the dike face
                        ys_new_r = ys_hs + (sramp_new_r/betas) - (sramp_new_r/betas_hs)
                        dune_x['ramp_crest'] = ys_new_r - (sramp_new_r/betas)
                        dune_x['ramp_dike_temp'] = dune_x['seaward_crest']
                        dune_z['ramp_dike_temp'] = sramp_new_r
                        dune_x['ramp_toe'] = ys_new_r
                        dune_z['ramp_crest'] = sramp_new_r
                        
        
                    else: #needs new calculationgs
                        vol_to_fill_top = ((ysprim_hs - dune_x['seaward_crest'])**2)*betas*0.5 #volume to fill the top of the dike and stop the intersection
                        vol_to_fill_face = ((((s_hs)**2)/(2*betas)))-((ys_hs-ysprim_hs)*(s_hs)*0.5) #volume to fill the face of the dike
                        vol_extra = vrampnew - vol_to_fill_face #extra volume that needs to be accomodated above the dike
        
                        if vol_extra < vol_to_fill_top: #thus not enough to fill everything and stop connection with dike
        
                            def f_ramp_dike_top(sramp_potential):
                                toebound = ysprim_hs
                                crestbound = dune_x['seaward_crest']
                                sramp_prim_x = toebound - (sramp_potential/betas)
                                dune_x_crest = np.array([crestbound, sramp_prim_x, toebound])
                                dune_z_crest = np.array([sramp_potential - s_hs, sramp_potential - s_hs, 0])
        
                                # volume between crest line and substrate
                                V = np.trapz(dune_z_crest, dune_x_crest)
                        
                                return V - vol_extra
                            
                            # solve for sramp_potential that gives the target volume
                            sramp_new_r = brentq(f_ramp_dike_top, 0, dune_z['seaward_crest']-s_hs) + s_hs # find the new ramp height using the root finder padding only the top of the dike
                            dune_x['ramp_crest'] = ysprim_hs - (sramp_new_r/betas)
                            dune_z['ramp_crest'] = sramp_new_r
                            
                        else: #enough to fill everything and stop connection with dike
                            ysnew = sramp_new_r / betas + dune_x['seaward_crest']
                            dune_x['ramp_crest'] = dune_x['seaward_crest']
                            dune_z['ramp_crest'] = sramp_new_r + dfoot
                            dune_x['ramp_dike_temp'] = None
                            dune_z['ramp_dike_temp'] = None
                            dune_x['ramp_toe'] = ysnew
                            dune_z['ramp_toe'] = dfoot
        
                else:
                    sramp_new_r = s  # ramp full
                    vextra = vrampnew - vrampmax # extra volume that needs to be accomodated by moving the dune toe and ysprim
                    
                    ys_new_r = (vextra / s) + (ys + ((sramp_new_r - sramp) / betas)) # get the new seaward position of the ramp
                    ysprim_new_r = ysprim + (vextra / s)  # get the new crest position of the ramp 
                    
                    dune_x['ramp_toe'] = ys_new_r
                    dune_z['ramp_crest'] = sramp_new_r+dfoot
                    dune_x['ramp_crest'] = ysprim_new_r
                    dune_x['seaward_crest'] = ysprim_new_r
                    
        
            else: #continue with normal procedure without hard structure interaction
        
                if vrampnew < vrampmax: # check if the volume in the front of dune is below the maximum
                
                    def f_ramp(sramp_potential):
                        return ramp_volume(sramp_potential, dune_z['ramp_crest']-dfoot, dune_x['ramp_crest'], dune_x['ramp_toe'], betas, dfoot, substrate_x, substrate_z) - vrampnew
                
                    # solve for sramp_potential that gives the target volume
                    sramp_new_r = brentq(f_ramp, dune_z['ramp_crest']-dfoot, dune_z['seaward_crest']-dfoot) # find the new ramp height using the root finder
                    ys_new_r = dune_x['ramp_crest'] + (sramp_new_r / betas)
                    #ys_new_r = ys + ((sramp_new_r-sramp) / betas)
            
                else: # if the volume exceeds the maximum, then move the role thing with fixed ramp height and slope, but moving ysprim.
                    sramp_new_r = dune_z['seaward_crest']-dfoot  # ramp full
                    vextra = vrampnew - vrampmax # extra volume that needs to be accomodated by moving the dune toe and ysprim
                    
                    ys_new_r = (vextra / (dune_z['seaward_crest']-dfoot)) + (ys + ((sramp_new_r - sramp) / betas)) # get the new seaward position of the ramp
                    ysprim_new_r = ysprim + (vextra / s)  # get the new crest position of the ramp 
                
                # update dune array with new ramp and toe position
                if vrampnew < vrampmax:
                    dune_x['ramp_toe'] = ys_new_r
                    dune_x['ramp_toe_temp'] = ys_new_r
                    dune_z['ramp_crest'] = sramp_new_r+dfoot
                    dune_z['ramp_dike_temp'] = sramp_new_r+dfoot
                else:
                    dune_x['seaward_crest'] = ysprim_new_r
                    dune_x['ramp_crest'] = ysprim_new_r
                    dune_x['ramp_dike_temp'] = ysprim_new_r
                    dune_x['ramp_toe'] = ys_new_r
                    dune_x['ramp_toe_temp'] = ys_new_r
                    dune_z['ramp_crest'] = sramp_new_r+dfoot
                    dune_z['ramp_dike_temp'] = sramp_new_r+dfoot
                    #dune_z['crest_temp'] = sramp_new_r+dfoot
                    #dune_x['crest_temp'] = ysprim_new_r
                
                # calculate total initial volume above the substrate
        total_volume_ramp = total_volume([dune_x[k] for k in list(dune_x.keys())[1:8]], [dune_z[k] for k in list(dune_z.keys())[1:8]], substrate_x, substrate_z)
        
        # calculate total initial volume above the substrate
        if dune_x['ramp_dike_temp'] is not None:
            ramp_idx = ['ramp_dike_temp', 'ramp_crest', 'ramp_toe']
        else:
            ramp_idx = ['ramp_crest', 'ramp_toe']
        
        total_volume_ramp_ramp = total_volume([dune_x[k] for k in ramp_idx], [dune_z[k] for k in ramp_idx], substrate_x, substrate_z)
        total_volume_crest_ramp = total_volume([dune_x[k] for k in crest_idx], [dune_z[k] for k in crest_idx], substrate_x, substrate_z)
        total_volume_lee_ramp = total_volume([dune_x[k] for k in lee_idx], [dune_z[k] for k in lee_idx], substrate_x, substrate_z)
        
        ################## Crest ##################
        vlee_lateral = 0
        detach_tag = 0
        
        if qs_crest > 0:
            # change in the crest
            deltavcrest = qs_crest # how much sand is deposited on the crest
            
            # calculate the substrate under the crest
            vcrestnew = total_volume_crest_ramp + deltavcrest # new volume of the crest after accretion
            
            # check if needed to use temp of lee_crest due to detached lee_side
            
            if dune_x['lee_crest'] <= substrate_x[2] and dune_z['lee_crest'] <= substrate_z[2]:
                # calculate maximum volume up till the crest becomes triangular (when ysprim=ylprim)
                s_max_triangular = (betal * betas * (dune_x['seaward_crest'] - dune_x['crest_temp'])) / (betal + betas)
                vcrest_triangular = 0.5 * (dune_x['seaward_crest'] - dune_x['crest_temp']) * s_max_triangular #max volume for triangular shape
                       
            else:
                # calculate maximum volume up till the crest becomes triangular (when ysprim=ylprim)
                s_max_triangular = (betal * betas * (dune_x['seaward_crest'] - dune_x['lee_crest'])) / (betal + betas)
                vcrest_triangular = 0.5 * (dune_x['seaward_crest'] - dune_x['lee_crest']) * s_max_triangular #max volume for triangular shape
        
            if deltavcrest <= vcrest_triangular:
                # calculate the pontial new dune height, ysprim and ylprim based on geometry and relation with the substrate
                # Define function for volume as function of s potential (crest height)
                def f_crest_triangular(s_potential):
        
                    ystemp = dune_x['seaward_crest'] # fix it but transform it in the new "foot" so to account for the increase as a trapezoid, not block
                    yltemp = dune_x['lee_crest']
        
                    ysprim_potential = ystemp - ((s_potential-s)/betas)
                    ylprim_potential = yltemp + ((s_potential-s)/betal)
            
                    # linear crest line: from (ylprim_potential, s_potential) to (ysprim_potential, s_potential)
                    dune_x_crest = np.array([yltemp, ylprim_potential, ysprim_potential, ystemp])
                    dune_z_crest = np.array([s + dfoot, s_potential + dfoot, s_potential + dfoot, s + dfoot])
            
                    n_per_segment = 100  # how many points per segment
                    x_dense = np.concatenate([np.linspace(dune_x_crest[i], dune_x_crest[i+1], n_per_segment, endpoint=False) for i in range(len(dune_x_crest)-1)] + [[dune_x_crest[-1]]])
                    z_dense = np.concatenate([np.linspace(dune_z_crest[i], dune_z_crest[i+1], n_per_segment, endpoint=False) for i in range(len(dune_z_crest)-1)] + [[dune_z_crest[-1]]])
            
                    # substrate under crest
                    z_sub_crest = np.interp(x_dense, substrate_x, substrate_z)
            
                    # volume between crest line and substrate
                    V = trapezoid_positive(z_dense - z_sub_crest, x_dense)
                    
                    if s_potential == s:
                        return -vcrestnew
                    else:
                        return V - vcrestnew
                
                def f_crest_triangular_detach(s_potential):
        
                    ystemp = dune_x['seaward_crest'] # fix it but transform it in the new "foot" so to account for the increase as a trapezoid, not block
                    #yltemp = dune_x['crest_temp']
                    yltemp = substrate_x[2] # fixed at dike position
        
                    ysprim_potential = ystemp - ((s_potential-s)/betas)
                    ylprim_potential = yltemp + ((s_potential-s)/betal)
            
                    # linear crest line: from (ylprim_potential, s_potential) to (ysprim_potential, s_potential)
                    dune_x_crest = np.array([yltemp, ylprim_potential, ysprim_potential, ystemp])
                    dune_z_crest = np.array([s + dfoot, s_potential + dfoot, s_potential + dfoot, s + dfoot])
            
                    n_per_segment = 100  # how many points per segment
                    x_dense = np.concatenate([np.linspace(dune_x_crest[i], dune_x_crest[i+1], n_per_segment, endpoint=False) for i in range(len(dune_x_crest)-1)] + [[dune_x_crest[-1]]])
                    z_dense = np.concatenate([np.linspace(dune_z_crest[i], dune_z_crest[i+1], n_per_segment, endpoint=False) for i in range(len(dune_z_crest)-1)] + [[dune_z_crest[-1]]])
            
                    # substrate under crest
                    z_sub_crest = np.interp(x_dense, substrate_x, substrate_z)
            
                    # volume between crest line and substrate
                    V = trapezoid_positive(z_dense - z_sub_crest, x_dense)
            
                    if s_potential == s:
                        return -vcrestnew
                    else:
                        return V - vcrestnew
                
                # check detached lee side
            
                if dune_x['lee_crest'] <= substrate_x[2] and dune_z['lee_crest'] <= substrate_z[2]:
                    s_new_c = brentq(f_crest_triangular_detach, s_hs, s_max_triangular+s)  # adjust 10 if necessary
                
                    # Calculate corresponding pivot points
                    ysprim_new_c = dune_x['seaward_crest'] - (s_new_c - s) / betas
                    ylprim_new_c = dune_x['crest_temp'] + (s_new_c - s) / betal
                    yl_new_c = dune_x['lee_foot']
                    z_int_c = dune_z['lee_foot']
                    ys_new_c = dune_x['ramp_toe']
                    vlee_lateral = 0.0
                
                    #check volume
                    h = s_new_c - s
                    b_back = dune_x['seaward_crest'] - dune_x['lee_crest']
                    b_front = ysprim_new_c - ylprim_new_c
                    V_check = 0.5 * (b_back + b_front) * h
                    
                    detach_tag = 1
                
                else:
                # Solve for s_potential that gives the target volume
                    s_new_c = brentq(f_crest_triangular, s, s_max_triangular+s)  # adjust 10 if necessary
            
                    # Calculate corresponding pivot points
                    ysprim_new_c = dune_x['seaward_crest'] - (s_new_c - s) / betas
                    ylprim_new_c = dune_x['lee_crest'] + (s_new_c - s) / betal
                    yl_new_c = dune_x['lee_foot']
                    z_int_c = dune_z['lee_foot']
                    ys_new_c = dune_x['ramp_toe']
                    vlee_lateral = 0.0
                
                    #check volume
                    h = s_new_c - s
                    b_back = dune_x['seaward_crest'] - dune_x['lee_crest']
                    b_front = ysprim_new_c - ylprim_new_c
                    V_check = 0.5 * (b_back + b_front) * h
                    
                    detach_tag = 0
                
            else:
                # volume that cannot fit into trapezoid (triangular crest)
                vcrest_lateral = deltavcrest - vcrest_triangular
            
                # assume half goes to ramp side, half goes to lee side
                vramp_lateral = vcrest_lateral / 2
                vlee_lateral = vcrest_lateral / 2
            
                # height of filled triangular crest
                s_new_c = s + s_max_triangular
            
                # pivot points for trapezoid (before lateral redistribution) those two should be the same because it it a triangle
                ysprim_tri = dune_x['seaward_crest'] - (s_max_triangular / betas)
                ylprim_tri = dune_x['lee_crest'] + (s_max_triangular / betal)
        
        
                # --- lateral overflow mapping ---
                # triangular cross-section: area = 0.5 * base * height -> base = 2 * area / height
                delta_ys = vramp_lateral / s_new_c
            
                # update positions
                ys_new_c = dune_x['ramp_toe'] + delta_ys            # seaward dune foot
                ysprim_new_c = ysprim_tri + delta_ys      # seaward pivot
                
                if s_new_c < s_hs: #check if there is space to grow behind the dune
                    # check if there is space in the lee side
                    lee_to_fill = lee_volume(dune_x['lee_crest'], dune_x, dune_z, betal, substrate_x, substrate_z)
                
                else:
                    yl_new_c = dune_x['lee_foot']                       # lee foot remains
                    ylprim_new_c = ylprim_tri                  # do not update yet because it will be updated in the next block
            
                # bottom of dune remains at dfoot
                z_int_c = dune_z['lee_foot']
            
            # Update dune profile
            
            if detach_tag == 1: 
                dune_x['crest_temp'] = ylprim_new_c
                dune_z['crest_temp'] = s_new_c + dfoot
            else:
                dune_x['lee_foot'] = yl_new_c
                dune_z['lee_foot'] = z_int_c
                dune_x['lee_crest'] = ylprim_new_c
                dune_z['lee_crest'] = s_new_c + dfoot
                dune_x['crest_temp'] = ylprim_new_c
                dune_z['crest_temp'] = s_new_c + dfoot
                
            dune_x['seaward_crest'] = ysprim_new_c
            dune_z['seaward_crest'] = s_new_c + dfoot
            dune_z['ramp_dike_temp'] = s_new_c + dfoot
            dune_x['ramp_dike_temp'] = ysprim_new_c
            dune_x['ramp_crest'] = ysprim_new_c
            dune_z['ramp_crest'] = s_new_c + dfoot
            dune_x['ramp_toe'] = ys_new_c
            dune_z['ramp_toe'] = dfoot
            dune_x['ramp_toe_temp'] = ys_new_c #TEST
            dune_z['ramp_toe_temp'] = dfoot #TEST

        
        ################## update lee ##################
        if qs_lee + vlee_lateral > 0:
            
            deltavlee = qs_lee + vlee_lateral # how much sand is deposited on the lee side
            deltavlee_temp = None
            
            vlee_c      = lee_volume(dune_x['lee_crest'], dune_x, dune_z, betal, substrate_x, substrate_z)
            vlee_target = vlee_c + deltavlee
            
            if s_new_c < s_hs and dune_x['lee_crest'] > substrate_x[3]: #check if there is space to grow behind the dune
            
                # find yl assuming no dike
                yl_bare = dune_x['lee_crest'] - s_new_c / betal
                                    
                # use shapely to find intersection between new lee line and potential hard structure
                substrate_line = LineString(zip(substrate_x, substrate_z))
                # lee line: from lee crest (ylprim, s_new) to theoretical flat toe (yl_flat, dfoot)
                lee_line = LineString([(dune_x['lee_crest'], s_new_c+dfoot), (yl_bare, dfoot)])
                # crest line: cres intersection
                crest_line = LineString([(dune_x['lee_crest'], s_new_c+dfoot), (ysprim_hs, s_new_c+dfoot)]) #only until ysprim because is the potential crossing when there is space behind
                # find intersection(s)
                intersection = substrate_line.intersection(lee_line)
                intersection_crest = substrate_line.intersection(crest_line)
        
                if not intersection.is_empty:
                    if intersection.geom_type == "Point":
                        yl_intersect, z_intersect = intersection.x, intersection.y
                    else:  # multiple intersections, take the shallowest one (max x)
                        points = [pt for pt in intersection.geoms if pt.geom_type == "Point"]
                        yl_intersect, z_intersect = max(points, key=lambda p: p.x).x, max(points, key=lambda p: p.x).y
               
                if not intersection_crest.is_empty:
                    if intersection_crest.geom_type == "Point":
                        yl_intersectcrest, z_intersectcrest = intersection_crest.x, intersection_crest.y
                    else:  # multiple intersections, take the shallowest one (max x)
                        points = [pt for pt in intersection_crest.geoms if pt.geom_type == "Point"]
                        yl_intersectcrest, z_intersectcrest = max(points, key=lambda p: p.x).x, max(points, key=lambda p: p.x).y
                        
                        
        
                # calculate the volume under the curve between ylprim and intersection point
                if dune_x['lee_crest'] > (ys_hs - (s_new_c+dfoot)/betas_hs): #check if there is space here
                    vol_to_fill_lee = (z_intersectcrest-z_intersect)*(dune_x['lee_crest']-yl_intersectcrest)*0.5
                else:
                    vol_to_fill_lee = 0
                    yl_intersectcrest = ys_hs - (s_new_c+dfoot)/betas_hs
                    z_intersectcrest = s_new_c+dfoot
                    
                if vol_to_fill_lee > deltavlee: # if the volume deltavlee is smaller than this, proceed normally
                
                    def f_lee(ylprim_potential):
                        return lee_volume(ylprim_potential, dune_x, dune_z, betal, substrate_x, substrate_z) - vlee_target
                    
                    # solve for ylprim_potential that gives the target volume
                    ylprim_new_l = brentq(f_lee, dune_x['lee_crest'] - 200, dune_x['lee_crest']) # find the new lee height using the root finder
                    yl_flat = ylprim_new_l - (dune_z['seaward_crest']-dfoot) / betal # get the new lee foot position assuming no hard structure
                    
                    # use shapely to find intersection between new lee line and potential hard structure
                    substrate_line = LineString(zip(substrate_x, substrate_z))
                    # lee line: from lee crest (ylprim, s_new) to theoretical flat toe (yl_flat, dfoot)
                    lee_line = LineString([(ylprim_new_l, dune_z['seaward_crest']), (yl_flat, dfoot)])
                    
                    # find intersection(s)
                    intersection = substrate_line.intersection(lee_line)
                    
                    if not intersection.is_empty:
                        if intersection.geom_type == "Point":
                            yl_new_l, z_int_l = intersection.x, intersection.y
                        else:  # multiple intersections, take the shallowest one (max x)
                            points = [pt for pt in intersection.geoms if pt.geom_type == "Point"]
                            yl_new_l, z_int_l = max(points, key=lambda p: p.x).x, max(points, key=lambda p: p.x).y
                    
                    else:
                        # no intersection with hard substrate → use flat geometry
                        yl_new_l = yl_flat
                        z_int_l = dfoot
                    
                    # Update dune profile
                    dune_x['lee_foot'] = yl_new_l
                    dune_z['lee_foot'] = z_int_l
                    dune_x['lee_crest'] = ylprim_new_l
                    dune_z['lee_crest'] = s_new_c + dfoot
                    
                    if dune_z['lee_foot'] > dune_z['back_ref']:
                        dune_x['back_ref'] = yl_new_l #change the reference to the dunefoot to avoid calculating extra volume from the link between ref and dunefoot above reference
                        dune_z['back_ref'] = z_int_l
                    
                else: # otherwise, fill this space and spread the remaining volume in the crest only
                
                    #adjust the points of intersection here for the volume calculations.
                    dune_x['lee_crest'] = yl_intersectcrest #adjust ylprim
                    dune_z['lee_crest'] = z_intersectcrest

                    dune_x['crest_temp'] = yl_intersectcrest #adjust ylprim
                    dune_z['crest_temp'] = z_intersectcrest
                            
                    dune_x['lee_foot'] = yl_intersectcrest
                    dune_z['lee_foot'] = dfoot
                
                    # calculate the substrate under the crest
                    leecrest_idx = ['lee_foot', 'lee_crest', 'seaward_crest']
                    vcrestlee = total_volume([dune_x[k] for k in leecrest_idx], [dune_z[k] for k in leecrest_idx], substrate_x, substrate_z)
                    vcrestleenew = vcrestlee + deltavlee - vol_to_fill_lee# new volume of the crest after accretion
                    
                    vtocrest = deltavlee - vol_to_fill_lee #volume to be spread over the crest
                    ylprim_int_i = yl_intersectcrest #point of contact of the crest after filling
        
                    def f_crest_nolee(ylprim_int_potential):
                        ystemp = dune_x['seaward_crest'] # fix it but transform it in the new "foot" so to account for the increase as a trapezoid, not block
                        s_potential = s_new_c + betas_hs*(yl_intersectcrest-ylprim_int_potential)
                        
                        ysprim_potential = ystemp - ((s_potential-s_new_c)/betas)
                
                        # linear crest line: from (ylprim_potential, s_potential) to (ysprim_potential, s_potential)
                        dune_x_crest = np.array([ylprim_int_potential, ysprim_potential, ystemp])
                        dune_z_crest = np.array([s_potential + dfoot, s_potential + dfoot, s_new_c + dfoot])
                
                        n_per_segment = 100  # how many points per segment
                        x_dense = np.concatenate([np.linspace(dune_x_crest[i], dune_x_crest[i+1], n_per_segment, endpoint=False) for i in range(len(dune_x_crest)-1)] + [[dune_x_crest[-1]]])
                        z_dense = np.concatenate([np.linspace(dune_z_crest[i], dune_z_crest[i+1], n_per_segment, endpoint=False) for i in range(len(dune_z_crest)-1)] + [[dune_z_crest[-1]]])
                
                        # substrate under crest
                        z_sub_crest = np.interp(x_dense, substrate_x, substrate_z)
                
                        # volume between crest line and substrate
                        V = trapezoid_positive(z_dense - z_sub_crest, x_dense)
                
                        return V - vcrestleenew 
                    
                    
                    if f_crest_nolee(ysprim_hs) > 0: # check if the volume overflows the dike or not or not
                        ylprim_int_new = brentq(f_crest_nolee, ysprim_hs, yl_intersectcrest) # find the new lee height using the root finder
                        s_new_lee = s_new_c + betas_hs*(yl_intersectcrest-ylprim_int_new)
                        ysprim_new_c = dune_x['seaward_crest'] - ((s_new_lee-s_new_c)/betas)
                        
                        #update the profile 
                        dune_x['back_ref'] = yl_new_c
                        dune_z['back_ref'] = s_hs+dfoot
                        dune_x['lee_foot'] = ylprim_int_new
                        dune_z['lee_foot'] = s_new_lee + dfoot
                        # dune_x['back_ref'] = yl_new_c
                        # dune_z['back_ref'] = s_hs+dfoot
                        # dune_x['lee_foot'] = yl_new_c
                        # dune_z['lee_foot'] = s_hs + dfoot
                        dune_x['lee_crest'] = ylprim_int_new
                        dune_z['lee_crest'] = s_new_lee + dfoot
                        dune_z['crest_temp'] = s_new_lee + dfoot
                        dune_x['crest_temp'] = ylprim_int_new
                        dune_x['seaward_crest'] = ysprim_new_c
                        dune_z['seaward_crest'] = s_new_lee + dfoot
                        dune_x['ramp_crest'] = ysprim_new_c
                        dune_z['ramp_crest'] = s_new_lee + dfoot
                        
                    else:
                        v_overflow = -f_crest_nolee(ysprim_hs) #volume to be redistributed on top of the dune and dike
                        
                        #use the code for when volume is bigger than the dike
                        
                        def f_crest_overdike_overflow(s_potential):
                        
                            ystemp = dune_x['seaward_crest']
                            yltemp = ylprim_hs
                        
                            ysprim_potential = ysprim_new_c - ((s_potential - s_hs) / betas)
                            ylprim_potential = yltemp + ((s_potential - s_hs) / betal)
                        
                            dune_x_crest = np.array([
                                yltemp,
                                ylprim_potential,
                                ysprim_potential,
                                ystemp
                            ])
                        
                            dune_z_crest = np.array([
                                s_hs + dfoot,
                                s_potential + dfoot,
                                s_potential + dfoot,
                                s_hs + dfoot
                            ])
                        
                            n_per_segment = 100
                            x_dense = np.concatenate([
                                np.linspace(dune_x_crest[i], dune_x_crest[i+1],
                                            n_per_segment, endpoint=False)
                                for i in range(len(dune_x_crest) - 1)
                            ] + [[dune_x_crest[-1]]])
                        
                            z_dense = np.concatenate([
                                np.linspace(dune_z_crest[i], dune_z_crest[i+1],
                                            n_per_segment, endpoint=False)
                                for i in range(len(dune_z_crest) - 1)
                            ] + [[dune_z_crest[-1]]])
                        
                            #z_sub_crest = np.interp(x_dense, substrate_x, substrate_z)
                        
                            V = trapezoid_positive(z_dense - (s_hs+dfoot), x_dense) #use s_hs because it is assumes filling up to that position
                        
                            return V - v_overflow

                        # Solve for s_potential that gives the target volume
                        s_new_lee = brentq(f_crest_overdike_overflow, s_hs, s_hs+5)  # adjust 10 if necessary
                        
                        # Calculate corresponding pivot points
                        ysprim_new_c = dune_x['seaward_crest'] - (s_new_lee - s_new_c) / betas
                        ylprim_int_new = ylprim_hs + (s_new_lee - (s_hs)) / betal
                        yl_new_c = ylprim_hs
                        
                        # # Update dune profile

                        dune_x['back_ref'] = 0
                        dune_z['back_ref'] = substrate_z[0]
                        dune_x['lee_foot'] = substrate_x[1]
                        dune_z['lee_foot'] = substrate_z[1] #as it filled everything, change the location to the back base of the dike.
                        dune_x['lee_crest'] = substrate_x[2]
                        dune_z['lee_crest'] = substrate_z[2]
                        dune_z['crest_temp'] = s_new_lee + dfoot
                        dune_x['crest_temp'] = ylprim_int_new
                        dune_x['seaward_crest'] = ysprim_new_c
                        dune_z['seaward_crest'] = s_new_lee + dfoot
                        dune_x['ramp_dike_temp'] = ysprim_new_c
                        dune_z['ramp_dike_temp'] = s_new_lee + dfoot
                        dune_x['ramp_crest'] = ysprim_new_c
                        dune_z['ramp_crest'] = s_new_lee + dfoot
                        
                    # # Update dune profile
                    # dune_x['back_ref'] = yl_new_c
                    # dune_z['back_ref'] = s_hs
                    # dune_x['lee_foot'] = ylprim_int_new
                    # dune_z['lee_foot'] = s_new_lee + dfoot
                    # dune_x['lee_crest'] = ylprim_int_new
                    # dune_z['lee_crest'] = s_new_lee + dfoot
                    # dune_z['crest_temp'] = s_new_lee + dfoot
                    # dune_x['crest_temp'] = ysprim_new_c
                    # dune_x['seaward_crest'] = ysprim_new_c
                    # dune_z['seaward_crest'] = s_new_lee + dfoot
                    # dune_x['ramp_crest'] = ysprim_new_c
                    # dune_z['ramp_crest'] = s_new_lee + dfoot
                        
                        
                    # if s_new_lee+dfoot > s_hs+dfoot: #check if the volume is bigger than the dike, then redistribuite the overflow.
                    
                    #     #v_overflow = trapezoid_positive(np.array([s_new_lee+dfoot, s_new_lee+dfoot, dfoot])-s_hs+dfoot, [ylprim_int_new, ysprim_new_c, dune_x['ramp_toe']]) #calculate the volume above the dike crest
                    #     def f_crest_overdike(s_potential):

                    #         ystemp = dune_x['seaward_crest'] # fix it but transform it in the new "foot" so to account for the increase as a trapezoid, not block
                    #         yltemp = ylprim_hs
                            
                    #         ysprim_potential = ysprim_new_c - ((s_potential-(s_new_c+dfoot))/betas)
                    #         ylprim_potential = yltemp + ((s_potential-(s_new_c+dfoot))/betal)
                    
                    #         # linear crest line: from (ylprim_potential, s_potential) to (ysprim_potential, s_potential)
                    #         dune_x_crest = np.array([yltemp, ylprim_potential, ysprim_potential, ystemp])
                    #         dune_z_crest = np.array([s_hs, s_potential, s_potential, s_new_c+dfoot])
                    
                    #         n_per_segment = 100  # how many points per segment
                    #         x_dense = np.concatenate([np.linspace(dune_x_crest[i], dune_x_crest[i+1], n_per_segment, endpoint=False) for i in range(len(dune_x_crest)-1)] + [[dune_x_crest[-1]]])
                    #         z_dense = np.concatenate([np.linspace(dune_z_crest[i], dune_z_crest[i+1], n_per_segment, endpoint=False) for i in range(len(dune_z_crest)-1)] + [[dune_z_crest[-1]]])
                    
                    #         # substrate under crest
                    #         z_sub_crest = np.interp(x_dense, substrate_x, substrate_z)
                    
                    #         # volume between crest line and substrate
                    #         V = trapezoid_positive(z_dense - z_sub_crest, x_dense)
                    
                    #         return V - vcrestleenew
                    
                    #     # Solve for s_potential that gives the target volume
                    #     s_new_lee = brentq(f_crest_overdike, s_hs, s_hs+5) - dfoot  # adjust 10 if necessary
                        
                    #     # Calculate corresponding pivot points
                    #     ysprim_new_c = dune_x['seaward_crest'] - (s_new_lee - s_new_c) / betas
                    #     ylprim_int_new = ylprim_hs + (s_new_lee - (s_hs-dfoot)) / betal
                    #     yl_new_c = ylprim_hs
                        
                    # Update dune profile


            else:
                
                # if dune moving on top of the dike, find the volume with the current height until reaching the back of the dike
                if dune_x['lee_foot'] > substrate_x[2] and dune_x['lee_foot'] < substrate_x[3]: # top of the dike
                    v_to_fill = (dune_x['lee_foot'] - substrate_x[2])*(dune_z['lee_crest']-substrate_z[2])
                    
                    #check if amount of volume is larger than needed:
                    if v_to_fill <= deltavlee:
                        #reduce the amount and use the overflow moving forward                       
                        deltavlee_temp = deltavlee.copy() #copy to change later to make sure we save the correct amount for conservation check
                        deltavlee = deltavlee - v_to_fill
                        
                        # update the parameters
                        dune_x['lee_foot'] = substrate_x[1]
                        dune_z['lee_foot'] = substrate_z[1]
                        dune_x['lee_crest'] = substrate_x[2]
                        dune_z['lee_crest'] = substrate_z[2]
                        dune_x['crest_temp'] = dune_x['lee_crest'] + ((dune_z['crest_temp']-dune_z['lee_crest'])/betal)
                        dune_x['back_ref'] = 0 #change the reference to the dunefoot to avoid calculating extra volume from the link between ref and dunefoot above reference
                        dune_z['back_ref'] = substrate_z[0]
                                 
                # Check the need to start a new lee portion behind the dike structure
                if (dune_x['lee_foot'] <= substrate_x[1] and dune_z['lee_foot'] == substrate_z[1]) and dune_z['lee_crest']<=substrate_z[2]:
                    
                    # check if the volume input exceed the needed to reconnect the lee side
                    x_leefoot_max = substrate_x[2] - (substrate_z[2]-substrate_z[1])/betal
                    
                    v_top = (dune_x['crest_temp']-substrate_x[2])*(dune_z['crest_temp']-substrate_z[2])/2
                    v_bottom_full = (substrate_x[2]-x_leefoot_max)*(substrate_z[2]-substrate_z[1])/2
                    v_dike = (substrate_x[2]-substrate_x[1])*(substrate_z[2]-substrate_z[1])/2
                    v_bottom = ((dune_x['lee_crest']-dune_x['lee_foot'])*(dune_z['lee_crest']-dune_z['lee_foot'])/2) - ((dune_x['lee_crest']-substrate_x[1])*(dune_z['lee_crest']-substrate_z[1])/2)
                    
                    v_to_fill = v_top + v_bottom_full - v_dike - v_bottom
                    
                    #check if the incoming sand fills the back of not
                    
                    if v_to_fill > deltavlee:
                    
                        #calculate the lee volume that is already there
                        
                        
                        def f_lee_behinddike(ylprim_potential):
                             
                            z_potential = (betal_hs*(ylprim_potential-substrate_x[1])) + substrate_z[1]
                            #z_potential = (betal_hs*(ylprim_potential-dune_x['lee_crest'])) + dune_z['lee_crest']
                            ylpotential = ylprim_potential - ((z_potential -dune_z['lee_foot']) / betal) 
    
                            # lee line
                            #dune_x_line = np.array([ylpotential, ylprim_potential, dune_x['crest_temp']]) #to calculate, keep the initial ylprim always otherwise u will have the same lee always!
                            #dune_z_line = np.array([dune_z['lee_foot'], z_potential, dune_z['crest_temp']]) # new dune height
    
                            dune_x_line = np.array([ylpotential, ylprim_potential, substrate_x[2]]) # use the pads from the substrate at the end otherwise will create fake volume due to the crest increase
                            dune_z_line = np.array([dune_z['lee_foot'], z_potential, substrate_z[2]]) # new dune height
    
                            n_per_segment = 100  # how many points per segment
                            x_dense = np.concatenate([np.linspace(dune_x_line[i], dune_x_line[i+1], n_per_segment, endpoint=False) for i in range(len(dune_x_line)-1)] + [[dune_x_line[-1]]])
                            z_dense = np.concatenate([np.linspace(dune_z_line[i], dune_z_line[i+1], n_per_segment, endpoint=False) for i in range(len(dune_z_line)-1)] + [[dune_z_line[-1]]])
                           
                            # substrate under lee side
                            z_sub_line = np.interp(x_dense, substrate_x, substrate_z)
    
                            # trapezoidal integration
                            return trapezoid_positive(z_dense - z_sub_line, x_dense) - deltavlee - v_bottom
                           
         
                        ylprim_new_l = brentq(f_lee_behinddike, dune_x['lee_foot'], substrate_x[2]) # find the new lee height using the root finder
                        z_ylprim = (betal_hs*(ylprim_new_l-substrate_x[1])) + substrate_z[1] #find the elevation at that point
                        yl_new_l = ylprim_new_l - ((z_ylprim - substrate_z[1])/betal)
                        
                        # Update dune profile
                        dune_x['lee_foot'] = yl_new_l
                        dune_z['lee_foot'] = substrate_z[1]
                        dune_x['lee_crest'] = ylprim_new_l
                        dune_z['lee_crest'] = z_ylprim
                    
                    else: #proceed without the influence of the dike anymore
                    
                        def f_lee(ylprim_potential):
                            return lee_volume(ylprim_potential, dune_x, dune_z, betal, substrate_x, substrate_z) - vlee_target
       
                        # solve for ylprim_potential that gives the target volume
                        ylprim_new_l = brentq(f_lee, dune_x['crest_temp'] - 200, dune_x['crest_temp']) # find the new lee height using the root finder
                        yl_flat = ylprim_new_l - (dune_z['crest_temp']-dune_z['lee_foot']) / betal # get the new lee foot position assuming no hard structure
      
                        # Update dune profile
                        dune_x['lee_foot'] = yl_flat
                        dune_x['lee_crest'] = ylprim_new_l
                        dune_z['lee_crest'] = dune_z['crest_temp'] #reatach the temp with the normal
                        dune_x['crest_temp'] = ylprim_new_l
                        
                        if dune_x['back_ref'] > dune_x['lee_foot']:
                            dune_x['back_ref'] = dune_x['lee_foot']
                            dune_z['back_ref'] = dune_z['lee_foot']
                        
                else:    
                
                    def f_lee(ylprim_potential):
                        return lee_volume(ylprim_potential, dune_x, dune_z, betal, substrate_x, substrate_z) - vlee_target
                    
                    # solve for ylprim_potential that gives the target volume
                    ylprim_new_l = brentq(f_lee, dune_x['lee_crest'] - 200, dune_x['lee_crest']) # find the new lee height using the root finder
                    yl_flat = ylprim_new_l - (dune_z['seaward_crest']-dune_z['lee_foot']) / betal # get the new lee foot position assuming no hard structure
                    
                    # use shapely to find intersection between new lee line and potential hard structure
                    substrate_line = LineString(zip(substrate_x, substrate_z))
                    # lee line: from lee crest (ylprim, s_new) to theoretical flat toe (yl_flat, dfoot)
                    lee_line = LineString([(ylprim_new_l, dune_z['seaward_crest']), (yl_flat, dune_z['lee_foot'])])
                    
                    # find intersection(s)
                    intersection = substrate_line.intersection(lee_line)
                    
                    if not intersection.is_empty:
                        if intersection.geom_type == "Point":
                            yl_new_l, z_int_l = intersection.x, intersection.y
                        else:  # multiple intersections, take the shallowest one (max x)
                            points = [pt for pt in intersection.geoms if pt.geom_type == "Point"]
                            yl_new_l, z_int_l = max(points, key=lambda p: p.x).x, max(points, key=lambda p: p.x).y
                    
                    else:
                        # no intersection with hard substrate → use flat geometry
                        yl_new_l = yl_flat
                        if dune_z['lee_foot'] == substrate_z[3]: #check if it is moving on top of the dike
                            z_int_l = substrate_z[1]
                        else:
                            z_int_l = dune_z['lee_foot']
                    
                    
                    ### ADD CONDITION WHERE WHEN CREST IS FILLED AND THE LEE IS PADDED TO THE DIKE, START FILLING THE LEE SIDE FROM BEHIND THE DIKE 
                    # Update dune profile
                    dune_x['lee_foot'] = yl_new_l
                    dune_z['lee_foot'] = z_int_l
                    dune_x['lee_crest'] = ylprim_new_l
                    dune_z['lee_crest'] = s_new_c + dfoot
                    dune_x['crest_temp'] = ylprim_new_l
                    dune_z['crest_temp'] = s_new_c + dfoot
                    
                    if dune_x['back_ref'] > dune_x['lee_foot']:
                        dune_x['back_ref'] = yl_new_l
                        dune_z['back_ref'] = z_int_l

            if deltavlee_temp is not None:
                deltavlee = deltavlee_temp
                deltavlee_temp = None
                    
    return dune_x, dune_z, deltavramp, erod_overwash, deltavramp_ero, deltavall

def subaq(DT, VOLBAR, HO, T, WFALL, SLR, BARNOUR, CB, BACTIVE):
    # Empirical coefficients
    ACOEFF = 1.333
    BCOEFF = 1.0
    DCOEFF = 0.56e-6
    MCOEFF = -0.5

    # Dean parameter
    DEAN = HO / (WFALL * T)

    # Deepwater wave steepness
    LO = 1.56 * T**2
    STP = HO / LO

    # Equilibrium bar volume
    VOLBAREQ = CB * DEAN**ACOEFF * STP**BCOEFF * LO**2

    # Response coefficient for bar volume
    LAMBDA2 = DCOEFF * DEAN**MCOEFF

    if VOLBAREQ - VOLBAR > 0:
        LAMBDA = LAMBDA2
    else:
        LAMBDA = 0.30 * LAMBDA2  # Adjust as needed

    # Calculate bar volume change
    DVOL = (VOLBAREQ - VOLBAR) * (1.0 - np.exp(-LAMBDA * DT))

    # Transport to (from) bar
    QBAR = DVOL / DT

    # Transport onshore from shoreface nourishment
    NOURCOEFF = 0.1
    DVOLNOUR = BARNOUR * (1.0 - np.exp(-NOURCOEFF * LAMBDA * DT))
    QNOUR = DVOLNOUR / DT

    # Bruun rule transport due to sea level rise
    QSLR = SLR * BACTIVE

    return QBAR, QNOUR, QSLR, DVOL



def bar(volbar, qbar, dt, vnourbar, barnour, qnour):
    volbar_new = volbar + qbar * dt
    barnour_new = barnour + vnourbar - qnour * dt

    volbar = volbar_new
    barnour = barnour_new

    if volbar < 0:
        print('BAR VOLUME EXCEEDED')

    return volbar, barnour


def FALVEL(TEMP, D50):
    SRATIO = 2.65
    NYVAL = np.zeros(9)
    TETA = np.zeros(9)
    NYVAL[0] = 1.79
    NYVAL[1] = 1.51
    NYVAL[2] = 1.31
    NYVAL[3] = 1.14
    NYVAL[4] = 1.00
    NYVAL[5] = 0.894
    NYVAL[6] = 0.799
    NYVAL[7] = 0.724
    NYVAL[8] = 0.658
    NYVAL = NYVAL * 1E-6

    for I in range(9):
        TETA[I] = (I + 1) * 5.0

    if (TEMP < TETA[0]) or (TEMP > TETA[8]):
        print('TEMPERATURE FOR VISCOSITY CALCULATION IS OUT OF RANGE')
        print(f'{TEMP}°C')

    # Linear interpolation
    NY = np.interp(TEMP, TETA, NYVAL)
    B = (SRATIO - 1) * 9.81 * D50 ** 3 / NY ** 2

    if B < 39:
        WFALL = (SRATIO - 1) * 9.81 * D50 ** 2 / 18 / NY
    elif 39 <= B < 1E4:
        WFALL = ((SRATIO - 1) * 9.81) ** 0.7 * D50 ** 1.1 / 6 / NY ** 0.4
    else:
        WFALL = ((SRATIO - 1) * 9.81 * D50 / 0.91) ** 0.5

    return WFALL



def nour(i, iv, volat, volad, volab, volrb, maxnour):
    vnourd = 0.0
    vnourb = 0.0
    vnourbar = 0.0

    if i == int(volat[int(iv) - 1]):
        vnourd = volad[int(iv) - 1]
        vnourb = volab[int(iv) - 1]
        vnourbar = volrb[int(iv) - 1]
        print(f'beach nourishment {vnourb:.2f}')
        print(f'bar nourishment {vnourbar:.2f}')
        
        if iv < maxnour:
            iv += 1

    return vnourd, vnourb, vnourbar, iv

