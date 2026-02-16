"""
Debug script to understand what's happening in morpho_update when dune stops growing.
"""
import sys
sys.path.insert(0, 'src')

import numpy as np
from scipy.optimize import brentq

# Simulate the conditions when dune stops growing
dune_z = {
    'seaward_crest': 9.0,
    'ramp_crest': 8.5001,
    'ramp_toe': 3.0
}

dune_x = {
    'seaward_crest': 130.0,
    'ramp_crest': 189.0,
    'ramp_toe': 244.0
}

substrate_z = [3.0, 3.0, 8.5, 8.5, 3.0, 3.0]
substrate_x = [0, 20, 83, 189, 200, 400]

# Parameters
s_hs = 5.5
betas = 0.1
ysprim_hs = 189
ys_hs = 200
ysprim = 189

# Small deltavramp (one timestep worth)
deltavramp = 0.0047  # Typical qs_ramp value

print("Simulating one timestep of deposition with deltavramp =", deltavramp)
print("="*80)

# Calculate vol_extra
vol_to_fill_face = ((((s_hs)**2)/(2*betas)))-((ys_hs-ysprim_hs)*(s_hs)*0.5)
vrampnew = 120 + deltavramp  # Assume current vramp is 120
vol_extra = vrampnew - vol_to_fill_face

print(f"vol_to_fill_face = {vol_to_fill_face:.4f}")
print(f"vrampnew = {vrampnew:.4f}")
print(f"vol_extra = {vol_extra:.4f}")

# Calculate vol_to_fill_top
vol_to_fill_top = (ysprim_hs - (2*(dune_z['seaward_crest']-dune_z['ramp_crest'])/betas))*((dune_z['seaward_crest']-dune_z['ramp_crest'])/2)
print(f"vol_to_fill_top = {vol_to_fill_top:.4f}")

if vol_extra < vol_to_fill_top:
    print("\nEntering f_ramp_dike_top case (vol_extra < vol_to_fill_top)")
    
    # Define f_ramp_dike_top
    def f_ramp_dike_top(sramp_potential):
        toebound = dune_x['ramp_crest']
        crestbound = dune_x['seaward_crest']
        sramp_prim_x = toebound - (sramp_potential/betas)
        dune_x_crest = np.array([crestbound, sramp_prim_x, toebound])
        dune_z_crest = np.array([sramp_potential, sramp_potential, 0])
        V = np.trapezoid(dune_z_crest, dune_x_crest)
        return V - vol_extra
    
    # Check function values at bounds
    lower_bound = 0
    upper_bound = dune_z['seaward_crest'] - substrate_z[3]
    
    print(f"\nBrentq bounds: [{lower_bound}, {upper_bound}]")
    print(f"f(lower_bound) = {f_ramp_dike_top(lower_bound):.6f}")
    print(f"f(upper_bound) = {f_ramp_dike_top(upper_bound):.6f}")
    
    # Check if signs are opposite
    if f_ramp_dike_top(lower_bound) * f_ramp_dike_top(upper_bound) > 0:
        print("\n❌ ERROR: Function has same sign at both bounds!")
        print("   Brentq will fail or return a bound value.")
        print("   This is why the dune stops growing!")
    else:
        sramp_potential_solved = brentq(f_ramp_dike_top, lower_bound, upper_bound)
        sramp_new_r = sramp_potential_solved + substrate_z[3]
        print(f"\nSolved: sramp_potential = {sramp_potential_solved:.6f}")
        print(f"        sramp_new_r = {sramp_new_r:.6f}")
        print(f"        Previous ramp_crest = {dune_z['ramp_crest']:.6f}")
        print(f"        Height increase = {sramp_new_r - dune_z['ramp_crest']:.6f}")

