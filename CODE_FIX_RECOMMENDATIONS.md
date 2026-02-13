# Code Fix Recommendations - CS Hybrid Model

**Purpose**: Provide specific, actionable code change directions for each identified issue  
**Status**: Recommendations only - NO CODE CHANGES MADE  
**Date**: 2026-02-13

---

## Overview

This document provides concrete directions for fixing each of the 14 identified issues. Each recommendation includes:
- **What to change**: Specific code locations and logic
- **How to change it**: Pseudocode and approach
- **Why**: Physical justification
- **Testing**: How to validate the fix
- **Difficulty**: Estimated complexity (Easy/Medium/Hard)

---

## HIGH PRIORITY FIXES

### Issue #1: Beach Volume Fake Changes (main.py lines 483-525)

**Current Problem**: Beach volume calculation changes when `ramp_toe` moves, creating fake volume transfers.

#### Recommended Fix Approach

**Option A: Add Explicit Geometric Transfer Tracking (RECOMMENDED)**

```python
# Add BEFORE morpho_update (around line 483)
# Store old ramp_toe position
ramp_toe_old = dune_x['ramp_toe']

# Calculate old beach volume with old boundary
vbeach_old = ((dune_x['shoreline'] - ramp_toe_old) - p.bl) / p.al
db_old = vbeach_old / (dune_x['shoreline'] - ramp_toe_old)
vbeachtot_old = ramp_toe_old * (p.dclos + p.dfoot) + (dune_x['shoreline'] - ramp_toe_old) * (p.dclos + db_old)

# ... morpho_update happens here (line 491)

# AFTER morpho_update (around line 516)
# Get new ramp_toe position
ramp_toe_new = dune_x['ramp_toe']

# Calculate geometric volume transfer due to boundary shift
# This is volume that moves between compartments due to changing integration domain
if ramp_toe_new != ramp_toe_old:
    # Calculate volume in the strip between old and new ramp_toe
    delta_ramp_toe = ramp_toe_new - ramp_toe_old
    
    if delta_ramp_toe > 0:  # ramp_toe moved seaward
        # Beach "loses" the strip between old and new positions
        # But this is just bookkeeping - sediment is now counted in dune
        geometric_transfer_from_beach = delta_ramp_toe * (p.dclos + p.dfoot)
    else:  # ramp_toe moved landward  
        # Beach "gains" the strip
        geometric_transfer_from_beach = delta_ramp_toe * (p.dclos + p.dfoot)  # negative value
else:
    geometric_transfer_from_beach = 0.0

# Now calculate beach volume INCLUDING both physical and geometric changes
vbeachtotnew = vbeachtot_old + \
               (-qwnet - qslr + qnour - qbar) * p.dt + \  # Physical fluxes
               vnourb[i] + \                               # Nourishment
               deltavramp_ero + \                          # Eroded sediment from dune
               (-geometric_transfer_from_beach)            # Geometric correction (negative because this volume went to dune accounting)

# Calculate new shoreline with geometric correction included
dune_x['shoreline'] = dune_x['ramp_toe'] + (p.al * (vbeachtotnew - dune_x['ramp_toe'] * (p.dfoot + p.dclos)) + p.bl) / (1 + p.al * p.dclos)
```

**Physical Justification**: 
- Separates geometric bookkeeping from physical sediment transport
- The geometric transfer is NOT real sediment - it's volume that switches compartments due to boundary movement
- By tracking it explicitly, we maintain mass conservation while allowing flexible boundaries

**Testing**:
```python
# Zero-flux test
# Set qwnet=0, qbar=0, qslr=0, qnour=0, vnourb=0, deltavramp_ero=0
# Allow ramp_toe to move via morpho_update
# Check: vbeachtotnew - vbeachtot_old should equal -geometric_transfer_from_beach
# Check: Total system volume unchanged
```

**Difficulty**: Medium  
**Estimated time**: 2-3 hours  
**Risk**: Low - adds explicit tracking without changing core logic

---

**Option B: Use Fixed Reference Frame (More Complex)**

```python
# Define a FIXED reference point that never moves
x_ref_toe = p.ys_initial  # Store initial ramp toe position as reference

# Calculate beach volume using fixed reference, then subtract dune overlap
def calculate_beach_volume_fixed_ref(x_ref_toe, x_shoreline, dune_profile, p):
    """
    Calculate beach volume from fixed reference to shoreline,
    then subtract any dune volume that overlaps this domain.
    """
    # Total volume from reference to shoreline
    v_total_domain = (x_shoreline - x_ref_toe) * (p.dclos + p.dfoot/2)
    
    # Subtract dune volume that's seaward of actual ramp_toe
    dune_overlap = 0.0
    if dune_x['ramp_toe'] > x_ref_toe:
        # Dune extends into the reference domain
        # Calculate volume of dune between x_ref_toe and ramp_toe
        dune_overlap = calculate_dune_volume_in_range(x_ref_toe, dune_x['ramp_toe'], dune_profile)
    
    v_beach_actual = v_total_domain - dune_overlap
    return v_beach_actual

# Use this in main loop
vbeach = calculate_beach_volume_fixed_ref(x_ref_toe, dune_x['shoreline'], (dune_x, dune_z), p)
```

**Physical Justification**: 
- Beach is everything from fixed reference to shoreline that's NOT dune
- No geometric artifacts because reference never moves
- More physically clear: beach = total - dune

**Difficulty**: Hard  
**Estimated time**: 1-2 days  
**Risk**: Medium - requires significant refactoring

---

### Issue #2: Ramp Volume Geometric Inflation (functions.py lines 1279-1505)

**Current Problem**: Ramp volume increases quadratically with height change, creating fake volume.

#### Recommended Fix Approach

**Add Explicit Volume Constraint**

```python
# In morpho_update function, accretion case (around line 1279)

if qs_ramp > 0:
    deltavramp = qs_ramp  # Actual sediment deposited
    
    # Current ramp volume
    vramp_current = ramp_volume(dune_z['ramp_crest']-dfoot, dune_z['ramp_crest']-dfoot, 
                                dune_x['ramp_crest'], dune_x['ramp_toe'], betas, dfoot, 
                                substrate_x, substrate_z)
    
    # TARGET volume = current + deposited sediment
    vramp_target = vramp_current + deltavramp  # This is the ONLY volume we should have
    
    # Now find new profile that gives EXACTLY vramp_target
    # Current code calculates geometric volume which may be > vramp_target
    
    # OPTION 1: Solve for new profile constrained to exact volume
    def volume_constraint(params):
        """
        Find profile parameters (height, positions) that give exactly vramp_target volume.
        params = [sramp_new, ysprim_new] or similar
        Returns: calculated_volume - vramp_target (should be zero)
        """
        sramp_new, ysprim_new = params
        ys_new = ysprim_new + sramp_new / betas
        
        v_calc = ramp_volume(sramp_new, sramp_new, ysprim_new, ys_new, 
                            betas, dfoot, substrate_x, substrate_z)
        
        return v_calc - vramp_target
    
    # Solve with constraint
    from scipy.optimize import fsolve, least_squares
    
    # Initial guess based on current profile
    initial_params = [dune_z['ramp_crest']-dfoot, dune_x['ramp_crest']]
    
    # Solve for parameters that give exact target volume
    solution = least_squares(volume_constraint, initial_params, 
                           bounds=([0, x_lower], [s_upper, x_upper]))
    
    sramp_new, ysprim_new = solution.x
    ys_new = ysprim_new + sramp_new / betas
    
    # Update profile
    dune_x['ramp_crest'] = ysprim_new
    dune_z['ramp_crest'] = sramp_new + dfoot
    dune_x['ramp_toe'] = ys_new
    
    # VERIFY: Calculate resulting volume
    v_result = ramp_volume(sramp_new, sramp_new, ysprim_new, ys_new, 
                          betas, dfoot, substrate_x, substrate_z)
    
    if abs(v_result - vramp_target) > 1e-6:
        print(f"WARNING: Volume constraint not satisfied: target={vramp_target}, result={v_result}")
```

**Physical Justification**:
- Deposited volume should EXACTLY equal profile volume increase
- No geometric inflation allowed
- If sediment is deposited, profile grows just enough to accommodate it
- The profile shape (height vs extension) is chosen to maintain target volume

**Alternative: Distributed Deposition**

```python
# Instead of calculating new profile geometry, distribute sediment across existing ramp

# Calculate height increase that would give correct volume
# For a ramp with fixed toe position and slope:
# dV = dh × (current_width + dh/(2*betas))
# This is approximately: dV ≈ dh × current_width for small dh

current_width = dune_x['ramp_toe'] - dune_x['ramp_crest']
dh_approx = deltavramp / current_width

# Increase height uniformly
dune_z['ramp_crest'] += dh_approx

# Adjust toe position to maintain slope
dune_x['ramp_toe'] = dune_x['ramp_crest'] + (dune_z['ramp_crest'] - dfoot) / betas

# Verify volume increase
vramp_new = ramp_volume(...)
actual_increase = vramp_new - vramp_current

if actual_increase > deltavramp * 1.1:  # More than 10% error
    # Iterate to refine
    scaling_factor = deltavramp / actual_increase
    dune_z['ramp_crest'] = dune_z['ramp_crest_old'] + dh_approx * scaling_factor
    dune_x['ramp_toe'] = dune_x['ramp_crest'] + (dune_z['ramp_crest'] - dfoot) / betas
```

**Testing**:
```python
# Pure deposition test
# Deposit exactly 10 m³/m on ramp
# Check: vramp_new - vramp_old should equal 10.0 ± 1e-6
# Check: No other compartment volumes change
```

**Difficulty**: Medium-Hard  
**Estimated time**: 4-6 hours  
**Risk**: Medium - changes core profile evolution logic

---

### Issue #3: Wrong Slope Calculations (functions.py lines 876, 878)

**Current Problem**: Slope calculated as `z/dx` instead of `(z_top - z_bottom)/dx`

#### Recommended Fix (EASIEST FIX!)

```python
# In deposition_dune function (around line 876-878)

# CURRENT (WRONG):
# slope_ramp = dune_z['ramp_crest']/dx_ramp if dx_ramp != 0 else 0
# slope_crest = 0
# slope_lee = dune_z['seaward_crest']/dx_lee if dx_lee != 0 else 0

# CORRECTED:
slope_ramp = (dune_z['ramp_crest'] - dune_z['ramp_toe']) / dx_ramp if dx_ramp != 0 else 0
slope_crest = 0  # Crest is flat by definition, this is correct
slope_lee = (dune_z['lee_crest'] - dune_z['lee_foot']) / dx_lee if dx_lee != 0 else 0

# Note: These are now TRUE SLOPES (rise over run) instead of elevation/distance
```

**Physical Justification**:
- Slope = vertical change / horizontal change = (z_high - z_low) / dx
- Previous calculation used absolute elevation, which is wrong
- This affects wind acceleration/deceleration calculations
- Correct slope gives correct wind shear distribution

**Additional Change Needed - Slope Sign Convention**:

```python
# For lee slope, need to be careful about sign
# If lee slopes downward (typical), slope should be negative
# But we might want absolute value for calculations

# Check coordinates: does x increase landward or seaward?
# If seaward (lee_foot < lee_crest in x):
slope_lee = (dune_z['lee_crest'] - dune_z['lee_foot']) / dx_lee  # Positive if crest higher

# If landward (lee_foot > lee_crest in x):  
slope_lee = (dune_z['lee_crest'] - dune_z['lee_foot']) / dx_lee  # Negative if crest higher

# For wind calculations, we often want absolute value of lee slope
slope_lee_abs = abs(slope_lee)

# Update line 905 accordingly:
tau_lee = max(alpha_lee * (I_lee + beta_lee * -slope_lee_abs) + 1, 0.1)
# Note: The negative sign is already there, meaning lee side has deceleration
```

**Testing**:
```python
# Verify slopes are reasonable
# For typical dune:
# - Seaward slope: 5-15% (0.05-0.15)
# - Lee slope: 10-30% (0.10-0.30)
# 
# After fix, print and check:
print(f"Ramp slope: {slope_ramp:.3f} ({slope_ramp*100:.1f}%)")
print(f"Lee slope: {slope_lee:.3f} ({slope_lee*100:.1f}%)")
# 
# These should be physically reasonable
# Old values would be ~0.17 (too high) using elevation
# New values should be ~0.11 (correct) using height difference
```

**Difficulty**: Easy  
**Estimated time**: 30 minutes  
**Risk**: Very low - simple formula correction

---

## MEDIUM-HIGH PRIORITY FIXES

### Issue #4: Double Counting in Beach Volume Update (main.py lines 484-486, 517)

**Current Problem**: Eroded sediment may be counted twice - once geometrically and once explicitly.

#### Recommended Fix

```python
# The issue is that deltavramp_ero is added explicitly (line 517)
# But the boundary shift also implicitly includes this volume

# OPTION: Remove explicit addition, rely only on geometric transfer

# Around line 517, MODIFY:
# vbeachtotnew = vbeachtot_old + (-qwnet - qslr + qnour - qbar) * p.dt + vnourb[i] + deltavramp_ero

# Calculate ramp toe movement due to erosion
ramp_toe_old = dune_x['ramp_toe']  # Before morpho_update
# ... morpho_update happens ...
ramp_toe_new = dune_x['ramp_toe']  # After morpho_update

# If erosion moved ramp_toe landward, beach width increases automatically
# This geometric increase should account for the eroded volume
# So we might NOT need to add deltavramp_ero explicitly

# CHECK: Does morpho_update move ramp_toe in proportion to eroded volume?
# If yes, geometric transfer handles it
# If no, need explicit term

# RECOMMENDED: Add diagnostic to check
geometric_beach_width_increase = (ramp_toe_new - ramp_toe_old)  # Negative if moved landward
geometric_volume_to_beach = -geometric_beach_width_increase * (p.dclos + p.dfoot)

# Compare with explicit erosion
print(f"Geometric transfer: {geometric_volume_to_beach:.2f} m³/m")
print(f"Explicit erosion: {deltavramp_ero:.2f} m³/m")
print(f"Ratio: {geometric_volume_to_beach / max(deltavramp_ero, 1e-6):.2f}")

# If ratio ≈ 1.0, we're double counting!
# Solution: Only use geometric transfer OR explicit, not both

# CORRECTED CODE:
if abs(geometric_volume_to_beach - deltavramp_ero) < 0.1 * deltavramp_ero:
    # Geometric and explicit are similar - double counting risk
    # Use only geometric (don't add deltavramp_ero)
    vbeachtotnew = vbeachtot_old + (-qwnet - qslr + qnour - qbar) * p.dt + vnourb[i]
else:
    # They're different - need explicit term
    vbeachtotnew = vbeachtot_old + (-qwnet - qslr + qnour - qbar) * p.dt + vnourb[i] + deltavramp_ero
```

**Difficulty**: Medium  
**Estimated time**: 2-3 hours  
**Risk**: Medium - requires understanding of coupling between erosion and boundary movement

---

### Issue #5: Complex Erosion Logic Edge Cases (functions.py lines 1020-1260)

**Current Problem**: Many conditional branches with potential gaps.

#### Recommended Fix

```python
# Add comprehensive validation and edge case handling

def morpho_update_erosion_safe(dune_x, dune_z, substrate_x, substrate_z, qs, ql, ...):
    """
    Morphological update with erosion - safer version with validation.
    """
    
    # Input validation
    assert qs >= 0, "Front erosion should be positive"
    assert ql >= 0, "Lee erosion should be positive"
    
    erod_total = (qs + ql) * dt
    erod_ramp = qs * dt
    erod_overwash = ql * dt
    
    # Calculate available volume to erode
    available_volume = partial_volume(dune_x, dune_z, substrate_x, substrate_z,
                                     dune_x['lee_foot'], dune_x['ramp_toe'])
    
    # Limit erosion to available volume
    if erod_total > available_volume:
        print(f"WARNING: Erosion demand ({erod_total:.2f}) exceeds available volume ({available_volume:.2f})")
        # Scale down erosion proportionally
        scale_factor = available_volume / erod_total
        erod_ramp *= scale_factor
        erod_overwash *= scale_factor
        erod_total = available_volume
    
    # Check if dike is exposed (current logic around line 1044)
    if r + swl <= (s_hs + dfoot):  # Water below dike crest
        # Calculate where runup reaches on dike
        z_seg = substrate_z[3:5]
        x_seg = substrate_x[3:5]
        order = np.argsort(z_seg)
        runup_reach_dike_x = np.interp(r + swl, np.array(z_seg)[order], np.array(x_seg)[order])
        
        # Calculate where runup reaches on dune
        z_seg = [dune_z['ramp_toe'], dune_z['ramp_crest']]
        x_seg = [dune_x['ramp_toe'], dune_x['ramp_crest']]
        order = np.argsort(z_seg)
        runup_reach_dune_x = np.interp(r + swl, np.array(z_seg)[order], np.array(x_seg)[order])
        
        # DECISION LOGIC with explicit cases
        if runup_reach_dike_x >= runup_reach_dune_x:
            # Dike fully exposed - no dune erosion
            erod_ramp = 0
            print(f"Dike exposed at x={runup_reach_dike_x:.2f}, dune at x={runup_reach_dune_x:.2f}")
        elif runup_reach_dune_x <= dune_x['lee_foot']:
            # Runup reaches behind entire dune - maximum erosion
            erod_ramp = min(erod_ramp, available_volume)
            print(f"Runup reaches lee side - maximum erosion")
        else:
            # Partial erosion - calculate available volume
            volume_to_runup = partial_volume(dune_x, dune_z, substrate_x, substrate_z,
                                            runup_reach_dune_x, dune_x['ramp_toe'])
            erod_ramp = min(erod_ramp, volume_to_runup)
            print(f"Partial erosion up to x={runup_reach_dune_x:.2f}, volume={erod_ramp:.2f}")
    else:
        # Dike overtopped - different logic
        # ... (similar validation)
        pass
    
    # After calculating erod_ramp, find new profile
    if erod_ramp > 1e-6:  # Non-negligible erosion
        # Use root finding (current approach)
        try:
            def f_ramp_erosion(x):
                volume_potential = partial_volume(dune_x, dune_z, substrate_x, substrate_z, 
                                                 x, dune_x['ramp_toe'])
                return volume_potential - erod_ramp
            
            # Robust root finding with bounds checking
            x_lower = runup_reach_dike_x if 'runup_reach_dike_x' in locals() else dune_x['lee_foot']
            x_upper = dune_x['ramp_toe']
            
            # Check function values at bounds
            f_lower = f_ramp_erosion(x_lower)
            f_upper = f_ramp_erosion(x_upper)
            
            if f_lower * f_upper > 0:
                print(f"WARNING: Root not bracketed. f({x_lower})={f_lower:.3f}, f({x_upper})={f_upper:.3f}")
                # Use bisection or alternative
                if abs(f_lower) < abs(f_upper):
                    ysprim_new = x_lower
                else:
                    ysprim_new = x_upper
            else:
                ysprim_new = brentq(f_ramp_erosion, x_lower, x_upper)
                
        except ValueError as e:
            print(f"ERROR in root finding: {e}")
            # Fallback: estimate position
            ysprim_new = dune_x['ramp_toe'] - erod_ramp / (dune_z['ramp_crest'] - dune_z['ramp_toe'])
            print(f"Using fallback estimate: ysprim_new={ysprim_new:.2f}")
        
        # Update profile with validation
        # ... (existing update logic with bounds checking)
    
    # Return with diagnostic info
    return dune_x, dune_z, erod_ramp, erod_overwash, diagnostic_dict
```

**Difficulty**: Hard  
**Estimated time**: 1-2 days  
**Risk**: Medium - extensive refactoring of critical logic

---

### Issue #6: Ramp-Dike Interaction Geometry (functions.py lines 1291-1406)

**Current Problem**: Complex geometric assumptions that may not match reality.

#### Recommended Fix

```python
# Replace analytical geometry with numerical integration

# CURRENT (around line 1294):
# vrampmax_hs = ((((s_hs)**2)/(2*betas)))-((((s_hs)**2)/(2*betas_hs)))

# PROPOSED: Calculate actual fillable volume using profile integration

def calculate_fillable_volume_dike_front(dune_x, dune_z, substrate_x, substrate_z, 
                                         betas, s_hs, ysprim_hs, ys_hs, dfoot):
    """
    Calculate volume needed to fill front of dike from current profile to full connection.
    Uses actual geometry instead of analytical approximation.
    """
    
    # Current dune ramp extends to:
    current_toe = dune_x['ramp_toe']
    current_crest = dune_x['ramp_crest']
    current_height = dune_z['ramp_crest'] - dfoot
    
    # Target: fill to dike face
    # Create target profile that connects to dike
    
    # If filling to dike height s_hs:
    target_height = s_hs
    target_crest_x = ysprim_hs  # At dike crest position
    target_toe_x = target_crest_x + target_height / betas  # Based on dune slope
    
    # Create dense profile for current state
    n_points = 1000
    x_current = np.linspace(current_crest, current_toe, n_points)
    z_current = np.interp(x_current, [current_crest, current_toe], 
                         [dune_z['ramp_crest'], dune_z['ramp_toe']])
    
    # Create dense profile for target state  
    x_target = np.linspace(target_crest_x, target_toe_x, n_points)
    z_target = np.interp(x_target, [target_crest_x, target_toe_x],
                        [target_height + dfoot, dfoot])
    
    # Substrate under both profiles
    z_sub_current = np.interp(x_current, substrate_x, substrate_z)
    z_sub_target = np.interp(x_target, substrate_x, substrate_z)
    
    # Calculate volumes
    v_current = np.trapz(np.maximum(z_current - z_sub_current, 0), x_current)
    v_target = np.trapz(np.maximum(z_target - z_sub_target, 0), x_target)
    
    # Volume to fill
    v_fill = v_target - v_current
    
    return max(v_fill, 0.0)  # Can't be negative

# Use in morpho_update:
vrampmax_hs = calculate_fillable_volume_dike_front(dune_x, dune_z, substrate_x, substrate_z,
                                                   betas, s_hs, ysprim_hs, ys_hs, dfoot)
```

**Physical Justification**:
- Uses actual profile geometry instead of idealized triangles
- Accounts for substrate shape
- Works for any configuration
- More robust to edge cases

**Difficulty**: Medium-Hard  
**Estimated time**: 4-6 hours  
**Risk**: Medium - replaces analytical with numerical (slower but more robust)

---

## MEDIUM PRIORITY FIXES

### Issue #7: Fetch Distance with Moving Boundary (main.py lines 448-450)

**Current Problem**: Fetch changes when ramp_toe moves, creating artificial feedback.

#### Recommended Fix

```python
# Use absolute coordinates for fetch calculation

# CURRENT (line 448):
# yr = dune_x['ramp_toe'] + (1 - (r + swl) / p.dfoot) * (dune_x['shoreline'] - dune_x['ramp_toe'])
# bdry = yr - dune_x['ramp_toe']

# PROBLEM: yr calculation uses ramp_toe, so fetch changes when ramp_toe moves

# PROPOSED: Calculate fetch based on physical inundation extent

# Physical concept: Fetch is the distance over which wind can pick up sediment
# This depends on:
# 1. Where water reaches (wetted zone)
# 2. Where dry sand begins

# Water reaches to elevation (r + swl)
# Find x-coordinate where profile intersects this elevation

# Create profile from ramp_toe to shoreline (simplified as linear for now)
slope_beach = dune_z['ramp_toe'] / (dune_x['shoreline'] - dune_x['ramp_toe'])
elevation_shoreline = 0.0  # At mean sea level

# Find where water level intersects beach profile
if r + swl >= dune_z['ramp_toe']:
    # Water reaches beyond ramp toe - no dry beach
    bdry = 0.0
elif r + swl <= elevation_shoreline:
    # Water at shoreline or below - entire beach dry
    bdry = dune_x['shoreline'] - dune_x['ramp_toe']
else:
    # Water somewhere on beach
    # x_water = x where z = r + swl
    # For linear profile from (ramp_toe, z_ramp_toe) to (shoreline, 0):
    # z = z_ramp_toe - slope_beach * (x - ramp_toe)
    # Solve for x when z = r + swl:
    x_water = dune_x['ramp_toe'] + (dune_z['ramp_toe'] - (r + swl)) / slope_beach
    
    # Fetch is from water line to ramp toe
    bdry = x_water - dune_x['ramp_toe']
    
    # Ensure positive
    bdry = max(bdry, 0.0)

# Now fetch is based on actual inundation, not boundary position
yr = dune_x['ramp_toe'] + bdry  # For compatibility with existing code
f = bdry / abs(np.cos(np.radians(wdir)))

# Continue with existing aeolian transport calculation
qw = qwe * (1 - np.exp(-p.cwrate * f))
```

**Physical Justification**:
- Fetch depends on physical inundation extent, not arbitrary boundaries
- When ramp_toe moves but inundation doesn't change, fetch shouldn't change
- More realistic representation of sediment availability

**Difficulty**: Medium  
**Estimated time**: 2 hours  
**Risk**: Low - improves realism without breaking core logic

---

### Issue #9: Foreshore Slope Calculation (main.py line 350)

**Current Problem**: Uses elevation/width instead of true slope.

#### Recommended Fix

```python
# CURRENT (line 350):
# slope_foreshore = dune_z['ramp_toe']/(dune_x['shoreline']- dune_x['ramp_toe'])

# CORRECTED: Calculate slope between two points on actual beach face
# The foreshore is typically defined as the zone between low tide and high tide
# Or more generally, the active beach face

# Option 1: Use slope from shoreline to ramp toe
slope_foreshore = (dune_z['ramp_toe'] - 0.0) / (dune_x['ramp_toe'] - dune_x['shoreline'])
# Note: Assuming z=0 at shoreline (MSL), which should be defined

# Option 2: Use actual beach profile points if available
# If you have intermediate beach profile points:
# slope_foreshore = (z_upper_beach - z_lower_beach) / (x_upper - x_lower)

# Option 3: Make it clear what reference is being used
# Define explicitly:
z_reference = 0.0  # Mean sea level or other datum
x_reference = dune_x['shoreline']  # Where reference elevation is

slope_foreshore = (dune_z['ramp_toe'] - z_reference) / (dune_x['ramp_toe'] - x_reference)

# Add validation
if slope_foreshore < 0:
    print(f"WARNING: Negative foreshore slope ({slope_foreshore:.4f}). Check profile coordinates.")
    slope_foreshore = abs(slope_foreshore)

if slope_foreshore > 0.5:  # 50% slope is unrealistic for beach
    print(f"WARNING: Very steep foreshore slope ({slope_foreshore:.4f} = {slope_foreshore*100:.1f}%). Check profile.")

# The 0.158/0.03 factor in line 351 seems to be a calibration
# This should be reviewed based on the new slope definition
r = (0.158/0.03)*slope_foreshore*np.sqrt(h0 * LO)
```

**Difficulty**: Easy  
**Estimated time**: 1 hour  
**Risk**: Low - mainly clarification of existing calculation

---

### Issue #10: Shoreline Position Formula (main.py line 521)

**Current Problem**: No validation for extreme cases.

#### Recommended Fix

```python
# Add bounds checking and validation

# CURRENT (line 521):
# dune_x['shoreline'] = dune_x['ramp_toe'] + (p.al * (vbeachtotnew - dune_x['ramp_toe'] * (p.dfoot + p.dclos)) + p.bl) / (1 + p.al * p.dclos)

# IMPROVED with validation:

# Calculate new shoreline position
numerator = p.al * (vbeachtotnew - dune_x['ramp_toe'] * (p.dfoot + p.dclos)) + p.bl
denominator = 1 + p.al * p.dclos

shoreline_new = dune_x['ramp_toe'] + numerator / denominator

# Validate result
if shoreline_new <= dune_x['ramp_toe']:
    print(f"ERROR: Shoreline ({shoreline_new:.2f}) landward of ramp toe ({dune_x['ramp_toe']:.2f})")
    print(f"  vbeachtotnew = {vbeachtotnew:.2f}")
    print(f"  This indicates severe erosion beyond model assumptions")
    
    # Handle extreme erosion
    # Option 1: Set to minimum beach width
    min_beach_width = 10.0  # Minimum physically reasonable beach width (m)
    shoreline_new = dune_x['ramp_toe'] + min_beach_width
    print(f"  Setting shoreline to minimum width: {shoreline_new:.2f}")
    
    # Option 2: Allow shoreline to retreat to ramp toe (zero beach)
    # shoreline_new = dune_x['ramp_toe']
    
    # Update beach volume to match this constrained position
    vbeachtotnew = ((shoreline_new - dune_x['ramp_toe']) - p.bl) * (1/p.al)

# Check for unrealistic progradation
max_beach_width = 1000.0  # Maximum reasonable beach width (m)
if (shoreline_new - dune_x['ramp_toe']) > max_beach_width:
    print(f"WARNING: Extreme beach progradation ({shoreline_new - dune_x['ramp_toe']:.1f} m)")
    shoreline_new = dune_x['ramp_toe'] + max_beach_width

# Assign validated position
dune_x['shoreline'] = shoreline_new

# Calculate final beach volume for output (line 524)
vbeach = ((dune_x['shoreline'] - dune_x['ramp_toe'])-p.bl)/p.al

# Final check
if vbeach < 0:
    print(f"ERROR: Negative beach volume ({vbeach:.2f})")
    vbeach = 0.0
```

**Difficulty**: Easy  
**Estimated time**: 1 hour  
**Risk**: Low - adds safety checks

---

## LOW PRIORITY FIXES (Code Quality)

### Issue #11-14: Documentation and Validation

These are primarily about adding:
- Comments explaining sign conventions
- Input validation
- Bounds checking
- Better error messages

Example for Issue #14 (Intersection comment):

```python
# CURRENT comment (line 1244):
# # multiple intersections, take the shallowest one (max x)

# CORRECTED comment:
# Multiple intersections exist. Select the most landward one (maximum x-coordinate),
# which represents where the lee slope first contacts the substrate.
# This is typically at the toe of the dike structure.
```

---

## Implementation Strategy

### Phase 1: Quick Wins (Week 1)
1. ✅ Fix Issue #3 (slope calculations) - 30 min
2. ✅ Fix Issue #9 (foreshore slope) - 1 hour  
3. ✅ Fix Issue #10 (shoreline validation) - 1 hour
4. ✅ Add documentation fixes - 1 hour

**Total: ~4 hours, Low risk**

### Phase 2: Core Fixes (Week 2-3)
5. ✅ Fix Issue #1 (beach volume, Option A) - 3 hours
6. ✅ Fix Issue #7 (fetch distance) - 2 hours
7. ✅ Add validation in Issue #5 (erosion) - 1 day

**Total: ~2 days, Medium risk**

### Phase 3: Advanced Fixes (Week 4+)
8. ✅ Fix Issue #2 (ramp volume constraint) - 1 day
9. ✅ Fix Issue #6 (dike interaction) - 1 day
10. ✅ Comprehensive testing and validation - 2-3 days

**Total: ~1 week, Higher risk**

---

## Testing After Each Fix

### Required Tests
1. **Zero-flux test**: No transport, volume should be constant
2. **Pure deposition test**: Known input = known output
3. **Erosion test**: Conservation check
4. **Regression test**: Compare with pre-fix results on known cases

### Acceptance Criteria
- Mass conservation error < 0.1%
- No negative volumes
- Physically reasonable slopes and positions
- Results match previous version for simple cases
- Better results for extreme cases

---

## Additional Recommendations

### Add Diagnostic Mode

```python
# Add to main.py initialization
DIAGNOSTIC_MODE = True  # Set via config

if DIAGNOSTIC_MODE:
    # Create detailed log file
    diag_file = open('diagnostic_log.txt', 'w')
    
    # Log all volume changes
    def log_volume_change(step, compartment, old_vol, new_vol, reason):
        change = new_vol - old_vol
        diag_file.write(f"{step},{compartment},{old_vol:.6f},{new_vol:.6f},{change:.6f},{reason}\n")
    
    # Log all position changes
    def log_position_change(step, point, old_x, new_x, reason):
        change = new_x - old_x
        diag_file.write(f"{step},{point},{old_x:.6f},{new_x:.6f},{change:.6f},{reason}\n")

# Use throughout simulation
if DIAGNOSTIC_MODE:
    log_volume_change(i, 'beach', vbeach_old, vbeach, 'flux_update')
    log_position_change(i, 'ramp_toe', ramp_toe_old, ramp_toe_new, 'morpho_update')
```

### Add Assertions

```python
def assert_mass_conservation(compartments, fluxes, dt, tolerance=1e-6):
    """
    Check mass conservation for current timestep.
    """
    total_volume_change = sum(compartments['new']) - sum(compartments['old'])
    total_flux = sum(fluxes) * dt
    
    error = abs(total_volume_change - total_flux)
    
    if error > tolerance:
        raise ValueError(f"Mass conservation violated! Error = {error:.6e}")
```

---

## Summary of Recommendations

| Issue | Priority | Difficulty | Time | Risk | Fix Type |
|-------|----------|------------|------|------|----------|
| #1 | HIGH | Medium | 3h | Low | Add tracking |
| #2 | HIGH | Medium-Hard | 6h | Medium | Constraint |
| #3 | HIGH | Easy | 0.5h | Very Low | Formula fix |
| #4 | MED-HIGH | Medium | 3h | Medium | Logic check |
| #5 | MED-HIGH | Hard | 2d | Medium | Validation |
| #6 | MED-HIGH | Medium-Hard | 6h | Medium | Numerical |
| #7 | MEDIUM | Medium | 2h | Low | Realism |
| #8 | MEDIUM | Medium | 3h | Low | Numerical |
| #9 | MEDIUM | Easy | 1h | Low | Clarification |
| #10 | MEDIUM | Easy | 1h | Low | Validation |
| #11-14 | LOW | Easy | 2h | Very Low | Documentation |

**Total estimated time for all fixes**: 3-4 weeks

**Recommended minimum**: Fix #1, #2, #3 (HIGH priority) = ~1 week

---

**End of Code Fix Recommendations**

*Remember: Test after EACH fix, don't change multiple things at once!*
