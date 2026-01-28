'''
Docstring for examples.oosterover.run_example_oosteroever
This script runs the CS-model for the Oosteroever site with vegetation
and compares the results to data from Strypsteen et al. (2024).
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

####################### SETUP PATHS #########################

PROJECT_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../../dunefront-backend/models")
)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)
    
from cs_model.main import run_model

####################### RUN MODEL #########################

results = run_model("config_vegetation.yaml")

################### LOAD MODEL OUTPUT ####################

output_filename = 'cs_model_oosteroever.pkl'

with open(output_filename, 'rb') as f:
    results = pickle.load(f)

# Convert lists to numpy arrays for easier handling
out_strypsteen = {}

for key, value in results.items():
    out_strypsteen[key] = np.array(value)

## Load comparison data from Strypsteen et al. (2024)
with open("validation/vegetation_points_strypsteen2024.pkl", "rb") as f:
    vegetation_data = pickle.load(f)

veg_dates, veg_coverage = zip(*vegetation_data)
volume_trend = np.load("validation/volume_trend_strypsteen2024.npy")
vdune = out_strypsteen['vdune']
dune_x = out_strypsteen['dune_x']
dune_z = out_strypsteen['dune_z']
substrate_x = out_strypsteen['substrate_x']
substrate_z = out_strypsteen['substrate_z']
dfoot = 4.5  # dune foot elevation from Strypsteen et al. (2024)

########################## PLOT #############################

# Prepare data for plotting
time_model = out_strypsteen['time']   # time in seconds from model output
veg_model = veg_dates  # time of the measured vegetation coverage 

# Calculate the average of the vegetation modelled for plotting
spec1df_model = np.array(out_strypsteen['rhov_df'])
spec1dt_model = np.array(out_strypsteen['rhov_dt'])
spec1db_model = np.array(out_strypsteen['rhov_db'])

veg_model = 100*(spec1df_model + spec1dt_model + spec1db_model) / 3
veg_model = 100*(spec1dt_model)
    
# Find indices in model time corresponding to veg_dates
veg_dates_utc = pd.to_datetime(veg_dates)  # convert to datetime
start_date = pd.to_datetime("2021-01-11 00:00:00+00:00")
time_model_index = start_date + pd.to_timedelta(time_model, unit='s')

# Create 2-row subplot figure
fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

# Top subplot: Dune volume
axes[0].plot(time_model_index, vdune, color='green', alpha=0.7, linewidth=2, label='CS-model')
axes[0].plot(time_model_index, volume_trend, color='black', linewidth=2, label='Strypsteen et al. (2024)')
axes[0].set_ylabel("Dune Volume (m³/m)", fontsize=12)
axes[0].grid(True, linestyle='--', alpha=0.5)
axes[0].legend(fontsize=10, loc='best')

# Bottom subplot: Vegetation coverage
axes[1].plot(time_model_index, veg_model, color='green', alpha=0.7, linewidth=2)
axes[1].plot(veg_dates, veg_coverage, 'ko-', linewidth=2)
axes[1].set_ylabel("Vegetation Coverage (%)", fontsize=12)
axes[1].grid(True, linestyle='--', alpha=0.5)

# X-axis formatting (shared)
axes[1].xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1,5,9)))  # quarterly ticks
axes[1].xaxis.set_major_formatter(mdates.DateFormatter('%b-%y'))
plt.xticks(fontsize=10, rotation=0)
plt.yticks(fontsize=10)

plt.tight_layout()
plt.show()

################ SAVE FIGURE ####################
fig.savefig("oosteroever_vegetation_volume_comparison.png", dpi=300)