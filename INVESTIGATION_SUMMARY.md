# Investigation Summary: Ramp Crest Height Issue

## Problem Statement
The ramp_crest was not increasing properly when reaching the dike and entering the `f_ramp_dike_top` code block. With the user's parameters (only aeolian transport, no storms), the dune grew much slower than expected.

## Root Cause Identified

The issue is **NOT** primarily about missing `+ dfoot` terms (though some were missing). The core problem is in the **volume calculation logic** when the ramp reaches the dike top.

### What Happens:
1. **Before t=5000h**: Dune grows normally, ramp_crest increases from 9.0 m to 8.5 m
2. **At t~5000h**: Ramp_crest reaches 8.5001 m (dike top elevation = substrate_z[3] = 8.5 m)
3. **After t=5000h**: **Dune volume stops growing** despite continued aeolian deposition (~0.0047 m³/m per timestep)
4. **By t=8760h**: 127 m³/m of sediment deposited but NOT added to dune volume!

### Technical Details:

When `ramp_crest` is at the dike top, the code calculates:
- `vol_to_fill_face = 121.0 m³/m` (volume to fill dike face - a CONSTANT)
- `vrampnew = 120.0047 m³/m` (current ramp volume)
- `vol_extra = vrampnew - vol_to_fill_face = -0.9953 m³/m` (**NEGATIVE!**)

The negative `vol_extra` causes the `f_ramp_dike_top` function to have the same sign at both bounds:
- `f(0) = +0.9953` (positive)
- `f(0.5) = +29.2453` (positive)

This violates brentq's requirement for opposite signs, so it cannot find a valid solution. The geometry gets "stuck" and sediment continues to be deposited without updating the profile.

### Why vol_to_fill_face is Wrong:

`vol_to_fill_face` is calculated as if filling the dike face from scratch:
```python
vol_to_fill_face = ((((s_hs)**2)/(2*betas)))-((ys_hs-ysprim_hs)*(s_hs)*0.5)
```

But when the ramp is ALREADY at the dike top, this entire "face" is already filled! The calculation should account for the current ramp state, not assume it starts from zero.

## Fixes Applied (Partial):

1. **✓ Fixed lines 1431, 1433, 1455**: Added `+ dfoot` for proper height conversion
2. **✓ Fixed lines 1411, 1465**: Added missing `sramp_new_r` definitions in else blocks
3. **✓ Fixed numpy compatibility**: Updated `np.trapz` to handle both numpy 1.x and 2.x
4. **⚠️ Attempted fix for vol_extra issue**: Added `vol_extra > 0` condition, but this doesn't fully resolve the problem

## What Still Needs Fixing:

The `vol_to_fill_face` calculation (line 1374) needs to be refactored to account for the current ramp geometry state:

```python
# Current (WRONG when ramp is at dike top):
vol_to_fill_face = ((((s_hs)**2)/(2*betas)))-((ys_hs-ysprim_hs)*(s_hs)*0.5)

# Should be something like:
if dune_z['ramp_crest'] >= substrate_z[3]:  # Already at/above dike top
    vol_to_fill_face = 0  # No face left to fill
else:
    # Calculate only the UNFILLED portion
    remaining_height = substrate_z[3] - dune_z['ramp_crest']
    vol_to_fill_face = calculate_unfilled_face_volume(remaining_height, ...)
```

This requires understanding the full geometry and may affect other parts of the code.

## Test Results:

### Short test (720h):
- ✓ ramp_crest heights valid throughout
- ✓ Mass conservation: 0.0021% error (excellent)
- Dune experienced net erosion (expected for this scenario)

### Full test with user's parameters (8760h):
- ✓ ramp_crest heights valid until ~5000h
- ❌ Mass conservation: 65.98% error
- Dune volume: 192.48 → 313.48 m³/m (+121 m³/m)
- Expected: +248 m³/m
- **Missing: 127 m³/m of sediment!**

## Recommendation:

The `morpho_update` function needs a deeper refactoring of the dike-ramp interaction logic, particularly:
1. How `vol_to_fill_face` is calculated based on current geometry
2. The transition logic between different growth regimes
3. Proper handling when `vol_extra` becomes negative or very small

The current fixes improve code correctness but don't fully resolve the mass conservation issue when the ramp reaches the dike top with your specific parameters.
