'''
Docstring for examples.katwijk.run_example_katwijk
This script runs the CS-model for the Katwijk site with a hard structure
and compares the results to XBeach simulation results made by PhD Student Cem Sevindik (TU Delft)
'''

####################### IMPORTS #########################
import pandas as pd
import numpy as np
import pickle
from datetime import datetime, timedelta
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import sys
import os
import netCDF4 as nc

######################## SETUP PATHS #########################

PROJECT_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../../dunefront-backend/models")
)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from csmodel_hybrid.main import run_model

####################### RUN MODEL #########################

results = run_model("config_hardstruct.yaml", check_mass_conservation=True)

################### LOAD MODEL OUTPUT #################

output_filename = 'cs_model_katwijk.pkl'

with open(output_filename, 'rb') as f:
    results = pickle.load(f)

################## PROCESS DATA ####################

dune_x_i = results['dune_x'][0]
dune_z_i = results['dune_z'][0]
dune_x_10min = results['dune_x'][1]
dune_z_10min = results['dune_z'][1]
dune_x_20min = results['dune_x'][2]
dune_z_20min = results['dune_z'][2]
dune_x_60min = results['dune_x'][6]
dune_z_60min = results['dune_z'][6]
dune_x_120min = results['dune_x'][12]
dune_z_120min = results['dune_z'][12]
dune_x_240min = results['dune_x'][24]
dune_z_240min = results['dune_z'][24]
dune_x_360min = results['dune_x'][36]
dune_z_360min = results['dune_z'][36]
dune_x_480min = results['dune_x'][48]
dune_z_480min = results['dune_z'][48]
dune_x = results['dune_x'][-1]
dune_z = results['dune_z'][-1]
substrate_x = results['substrate_x']
substrate_z = results['substrate_z']

################## PLOTTING ####################

## load XBEACH simulation results for comparison
xbeach_file = "validation/xbeach/SWL_SB_With_SLR.nc"
ds = nc.Dataset(xbeach_file)

# load hard structure profile
hard_struct_df = pd.read_excel("validation/xbeach/bathymetry_Katwijk_xbeachCem.xlsx", sheet_name="hard_struct", header=None)
hard_struct_x = hard_struct_df.iloc[:, 0].values
hard_struct_z = hard_struct_df.iloc[:, 1].values

# make a plot where we have the profile evolution where I have the initial dune profile, 10 minutes, 20 minutes, 1 hour, 2 hour, 4 hours, 6 hour, 8 hour and final profile
x_var = ds.variables['globalx'][:]  # cross-shore distance
time_var = ds.variables['globaltime'][:]  # time in seconds
zb_var = ds.variables['zb'][:]  # bed elevation

x_data = x_var.data
time_data = time_var.data
zb_data = zb_var.data

# select time indices for the desired times (in seconds)
time_indices = [0, 600, 1200, 3600, 7200, 14400, 21600, 28800, time_data[-1]]
time_idx = [np.abs(time_data - t).argmin() for t in time_indices]
time_labels = ['Initial XBeach Profile', '10 min', '20 min', '1 hr', '2 hr', '4 hr', '6 hr', '8 hr', 'Final XBeach Profile']

#adjust x axis to match the reference plot (shift by 50 m)
x_data = (x_data.squeeze()-np.max(x_data.squeeze()))*-1
hard_struct_x = (hard_struct_x - np.max(x_data))*-1

# Define colors to match the reference plot
colors = ['blue', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'red']
linestyles = ['-', '-.', '-.', '-.', '-.', '-.', '-.', '-.', '-']
linewidths = [2, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 2]

fig, ax = plt.subplots(1, 2, figsize=(16, 6))

for idx in range(len(time_idx)):
    zb_profile = zb_data[time_idx[idx]] + 6.2  # sum 6.2 to get to the same datum as the CS
    ax[0].plot(x_data.squeeze(), zb_profile.squeeze(), label=time_labels[idx], 
             color=colors[idx], linestyle=linestyles[idx], linewidth=linewidths[idx])
# plot hard structure
ax[0].plot(hard_struct_x, hard_struct_z, label='Hard Structure', color='black', linestyle='--', linewidth=2.5)
xlim = (-22.910261266933176, 481.11548660559663)
ylim = (-0.41250000000000003, 8.6625)
ax[0].set_ylim(ylim)
ax[0].set_xlim(xlim)
ax[0].tick_params(axis='x', labelsize=14)
ax[0].tick_params(axis='y', labelsize=14)
ax[0].invert_xaxis()
ax[0].set_xlabel('Cross-shore distance (m)', fontsize=16)
ax[0].set_ylabel('Bed elevation (m)', fontsize=16)
ax[0].set_title('XBeach Simulation Profile Evolution', fontsize=16)
ax[0].legend(fontsize=12)
ax[0].grid(True)

# Reuse the same keys order for all profiles
keys = list(dune_x.keys())[1:10]

# Base profiles
ax[1].plot([dune_x_i[k] for k in keys], [dune_z_i[k] for k in keys], label='Initial CS Profile', color='blue', linewidth=2)

# In-between snapshots with distinct colors
profiles = [
    (dune_x_10min, dune_z_10min, '10 minutes'),
    (dune_x_20min, dune_z_20min, '20 minutes'),
    (dune_x_60min, dune_z_60min, '1 hour'),
    (dune_x_120min, dune_z_120min, '2 hours'),
    (dune_x_240min, dune_z_240min, '4 hours'),
    (dune_x_360min, dune_z_360min, '6 hours'),
    (dune_x_480min, dune_z_480min, '8 hours')
]
cmap = plt.get_cmap('tab10')
for i, (dx, dz, label) in enumerate(profiles):
    ax[1].plot([dx[k] for k in keys], [dz[k] for k in keys],
             label=label, color=cmap(i), linestyle='-.', linewidth=1.25)

# Final updated profile
ax[1].plot([dune_x[k] for k in keys], [dune_z[k] for k in keys], label='Final CS Profile', color='red', linewidth=2)
ax[1].plot(substrate_x[0:-1], substrate_z[0:-1], label='Hard Structure', color='black', linestyle='--', linewidth=2.5)

ax[1].set_xlabel('Cross-shore distance (m)', fontsize=16)
#ax[1].set_ylabel('Elevation (m)', fontsize=14)
ax[1].set_title('CS Model Profile Evolution', fontsize=16)
ax[1].set_yticklabels([])
ax[1].tick_params(axis='x', labelsize=14)
ax[1].legend(fontsize=12)
ax[1].grid(True)
ax[1].invert_xaxis()

plt.tight_layout()
plt.show()

################ SAVE FIGURE ####################
fig.savefig('katwijk_profile_evolution.png', dpi=300)