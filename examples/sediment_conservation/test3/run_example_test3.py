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
    os.path.join(os.path.dirname(__file__), "/Users/fgalifornisilv/Documents/GitHub/csmodel-hybrid/src/")
)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from cs_model.main import run_model

####################### RUN MODEL #########################

results = run_model("config_test3.yaml")

################### LOAD MODEL OUTPUT #################

output_filename = 'cs_model_test3.pkl'

with open(output_filename, 'rb') as f:
    results = pickle.load(f)

################## PROCESS DATA ####################

dune_x_i = results['dune_x'][0]
dune_z_i = results['dune_z'][0]
dune_x = results['dune_x'][-1]
dune_z = results['dune_z'][-1]
substrate_x = results['substrate_x']
substrate_z = results['substrate_z']


# ---- Transport time series ----
qs_transport = results['qs_crest'] + results['qs_lee'] + results['qs_ramp']
qbar = results['qbar']

# --- Helper functions ---
def extract_dune_section(dune_x, dune_z):
    """
    Returns arrays (x,z) only between ramp_toe and lee_foot.
    """

    # limits of dune body
    x_min = dune_x["lee_foot"]
    x_max = dune_x["ramp_toe"]

    # collect all pivots inside this interval
    xs, zs = [], []
    for key in dune_x.keys():
        x = dune_x[key]
        z = dune_z[key]

        if x_min <= x <= x_max:
            xs.append(x)
            zs.append(z)

    # sort by x
    xs = np.array(xs)
    zs = np.array(zs)

    order = np.argsort(xs)
    return xs[order], zs[order]


def profile_volume(dune_x, dune_z, substrate_x, substrate_z, dx=0.05):

    # common grid
    xmin = min(substrate_x)
    xmax = max(substrate_x)
    xg = np.arange(xmin, xmax+dx, dx)

    # interpolate dune and substrate
    zdune = np.interp(xg, dune_x, dune_z)
    zsub  = np.interp(xg, substrate_x, substrate_z)

    # sediment thickness above substrate
    hsed = np.maximum(zdune - zsub, 0)

    # integrate volume (m3/m)
    V = np.trapz(hsed, xg)

    return V

# --- Check dune volume conservation ---

x_dune0, z_dune0 = extract_dune_section(dune_x_i, dune_z_i)
x_dune1, z_dune1 = extract_dune_section(dune_x, dune_z)


V0 = profile_volume(
    x_dune0,
    z_dune0,
    substrate_x,
    substrate_z
)

V1 = profile_volume(
    x_dune1,
    z_dune1,
    substrate_x,
    substrate_z
)

print("Initial dune volume =", V0, "m3/m")
print("Final dune volume =", V1, "m3/m")

dV_geom = V1 - V0

print("Geometric dune volume change =", dV_geom, "m3/m")

V_deposited = np.sum(qs_transport)
V_eroded = np.sum(results['eroded_volume'])
dV_expected = V_deposited - V_eroded

residual = dV_geom - dV_expected
print("Residual (should be ~0):", residual)

relative_error = residual / dV_geom
print("Relative error:", relative_error)


################## PLOTTING ####################

plt.figure(figsize=(10, 4.8))

# -------------------------------
# Colors (colorblind-friendly)
# -------------------------------
c_init  = "#1f77b4"   # blue
c_final = "#d62728"   # red
c_sub   = "0.3"       # grey
c_fill_pos = "#2ca02c"  # green (deposition)
c_fill_neg = "#ff7f0e"  # orange (erosion)

# -------------------------------
# Extract pivot arrays
# -------------------------------
x_init = np.array(list(dune_x_i.values()))
z_init = np.array(list(dune_z_i.values()))

x_final = np.array(list(dune_x.values()))
z_final = np.array(list(dune_z.values()))

# Sort (important for interpolation/fill)
order_i = np.argsort(x_init)
order_f = np.argsort(x_final)

x_init, z_init = x_init[order_i], z_init[order_i]
x_final, z_final = x_final[order_f], z_final[order_f]

# -------------------------------
# Common interpolation grid
# -------------------------------
xg = np.linspace(
    max(x_init.min(), x_final.min()),
    min(x_init.max(), x_final.max()),
    500
)

z_init_g  = np.interp(xg, x_init, z_init)
z_final_g = np.interp(xg, x_final, z_final)

# -------------------------------
# Shaded difference (deposition/erosion)
# -------------------------------

# Deposition (final above initial)
plt.fill_between(
    xg, z_init_g, z_final_g,
    where=(z_final_g > z_init_g),
    interpolate=True,
    alpha=0.35,
    label="Deposition"
)

# Erosion (final below initial)
plt.fill_between(
    xg, z_init_g, z_final_g,
    where=(z_final_g < z_init_g),
    interpolate=True,
    alpha=0.35,
    label="Erosion"
)

# -------------------------------
# Plot profiles
# -------------------------------
plt.plot(x_init, z_init,
         color=c_init, lw=2.8,
         label="Initial profile")

plt.plot(x_final, z_final,
         color=c_final, lw=2.8,
         label="Final profile")

plt.plot(substrate_x[0:5],
         substrate_z[0:5],
         linestyle="--",
         color=c_sub,
         lw=2.0,
         label="Non-erodible substrate")

# -------------------------------
# Axis labels and title
# -------------------------------
plt.xlabel("Cross-shore distance (m)", fontsize=15)
plt.ylabel("Elevation (m)", fontsize=15)

plt.title("Cross-shore profile evolution and dune sediment conservation",
          fontsize=14)

# -------------------------------
# Tick label font size (publication)
# -------------------------------
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)

# -------------------------------
# Grid and legend
# -------------------------------
plt.grid(True, linestyle=":", linewidth=0.8, alpha=0.6)

plt.legend(loc="lower right",
           fontsize=11,
           frameon=True)

# -------------------------------
# Conservation annotation box
# -------------------------------
textstr = "\n".join([
    f"Initial dune volume: {V0:.1f} m³/m",
    f"Final dune volume:   {V1:.1f} m³/m",
    "",
    f"Residual:     {residual:.2e} m³/m",
    f"Rel. error:   {relative_error:.2%}"
])

plt.text(0.02, 0.05, textstr,
         transform=plt.gca().transAxes,
         fontsize=11.5,
         verticalalignment="bottom",
         bbox=dict(boxstyle="round,pad=0.45",
                   facecolor="white",
                   edgecolor="0.4",
                   alpha=0.95))

plt.tight_layout()

# -------------------------------
# Optional: Save high-res figure
# -------------------------------
# plt.savefig("dune_conservation.png", dpi=300)

plt.show()
