"""
Analyze the geometry created by f_ramp_dike_top to understand why it gets stuck.
"""
import numpy as np

# State when stuck
dune_z = {'seaward_crest': 9.0, 'ramp_crest': 8.5001}
dune_x = {'seaward_crest': 130.0, 'ramp_crest': 189.0}
substrate_z = [3.0, 3.0, 8.5, 8.5, 3.0, 3.0]
betas = 0.1

# Small vol_extra after reset
vol_extra = 0.0047

print("Analyzing f_ramp_dike_top geometry")
print("="*80)
print(f"Target vol_extra: {vol_extra:.6f} m³/m")
print(f"Current ramp_crest: z={dune_z['ramp_crest']:.4f} m, x={dune_x['ramp_crest']:.4f} m")
print(f"Seaward_crest: z={dune_z['seaward_crest']:.4f} m, x={dune_x['seaward_crest']:.4f} m")

# The f_ramp_dike_top function creates geometry:
toebound = dune_x['ramp_crest']  # 189.0
crestbound = dune_x['seaward_crest']  # 130.0

print(f"\nGeometry bounds:")
print(f"  toebound (ramp_crest x): {toebound:.2f} m")
print(f"  crestbound (seaward_crest x): {crestbound:.2f} m")
print(f"  x-span: {toebound - crestbound:.2f} m")

# Upper bound for brentq
upper_bound = dune_z['seaward_crest'] - substrate_z[3]  # 9.0 - 8.5 = 0.5
print(f"\nBrentq bounds: [0, {upper_bound:.4f}]")

# Test function at bounds
def volume_at_sramp(sramp_potential):
    sramp_prim_x = toebound - (sramp_potential/betas)
    dune_x_crest = np.array([crestbound, sramp_prim_x, toebound])
    dune_z_crest = np.array([sramp_potential, sramp_potential, 0])
    
    print(f"\n  sramp_potential = {sramp_potential:.4f}:")
    print(f"    sramp_prim_x = {toebound} - {sramp_potential}/{betas} = {sramp_prim_x:.2f}")
    print(f"    x_coords: [{dune_x_crest[0]:.2f}, {dune_x_crest[1]:.2f}, {dune_x_crest[2]:.2f}]")
    print(f"    z_coords: [{dune_z_crest[0]:.4f}, {dune_z_crest[1]:.4f}, {dune_z_crest[2]:.4f}]")
    
    # Check x-ordering
    if dune_x_crest[0] > dune_x_crest[1]:
        print(f"    ⚠️  WARNING: x[0] > x[1]! Trapezoid integration assumes increasing x!")
    if dune_x_crest[1] > dune_x_crest[2]:
        print(f"    ⚠️  WARNING: x[1] > x[2]! Trapezoid integration assumes increasing x!")
    
    V = np.trapezoid(dune_z_crest, dune_x_crest)
    f_val = V - vol_extra
    print(f"    Volume: {V:.6f} m³/m")
    print(f"    f(sramp) = V - vol_extra = {f_val:.6f}")
    
    return V, f_val

print("\nTesting at lower bound (0):")
V0, f0 = volume_at_sramp(0)

print("\nTesting at upper bound ({:.4f}):".format(upper_bound))
V_upper, f_upper = volume_at_sramp(upper_bound)

print("\n" + "="*80)
print("DIAGNOSIS:")
if f0 * f_upper > 0:
    print("❌ f(0) and f(upper) have the SAME SIGN!")
    print("   Brentq cannot find a root. It will likely return a bound value.")
    print("   This causes the height to not increase properly.")
else:
    print("✓ f(0) and f(upper) have OPPOSITE SIGNS")
    print("  Brentq should find a valid solution.")

print(f"\nExpected result: sramp_new_r = solution + substrate_z[3] = solution + {substrate_z[3]}")
print(f"This should be GREATER than current ramp_crest ({dune_z['ramp_crest']:.4f} m)")
print(f"for the dune to grow.")

