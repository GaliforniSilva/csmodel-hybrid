# Mass Conservation Analysis - CS Model

## Summary

✅ **Mass conservation issues have been identified and fixed**

The model was experiencing sediment mass balance errors of approximately 2.5% over a 24-hour simulation. After implementing the fixes described below, **mass is now conserved to within 0.0022% error**.

---

## Issues Found

### 1. Missing Bar-Beach Exchange Flux
**Location**: `src/cs_model/main.py`, line 507

**Problem**: The `qbar` flux (sediment exchange between bar and beach) was not included in the beach volume update equation.

**Impact**: Beach volume was not properly accounting for sediment deposited in or eroded from the offshore bar.

**Fix**: Added `-qbar` term to the beach volume equation:
```python
vbeachtotnew = vbeachtot_old + (-qwnet - qslr + qnour - qbar) * p.dt + vnourb[i] + deltavramp_ero
```

---

### 2. Geometric Discontinuity in Beach Volume Calculation
**Location**: `src/cs_model/main.py`, lines 470-503

**Problem**: Beach volume (`vbeachtot`) was calculated AFTER `morpho_update`, which updates the dune profile including the `ramp_toe` position. When `ramp_toe` moves due to erosion/accretion, the beach reference frame changes, creating a geometric discontinuity.

**Impact**: The ramp_toe movement implicitly changed beach geometry in a way that wasn't properly accounted for in the volume balance.

**Fix**: Calculate `vbeachtot_old` BEFORE calling `morpho_update`, using the OLD ramp_toe position:
```python
# Calculate beach volume BEFORE morpho_update (using old ramp_toe)
vbeach_old = ((dune_x['shoreline'] - dune_x['ramp_toe'])-p.bl)/p.al
db_old = vbeach_old / (dune_x['shoreline'] - dune_x['ramp_toe'])
vbeachtot_old = dune_x['ramp_toe']  * (p.dclos + p.dfoot) + (dune_x['shoreline'] - dune_x['ramp_toe']) * (p.dclos + db_old)

# Then call morpho_update...
# Then use vbeachtot_old for beach volume update
```

---

### 3. Eroded Sediment Accounting
**Location**: `src/cs_model/main.py`, line 507

**Problem**: When the dune erodes, the eroded sediment volume (`deltavramp_ero`) must be added to the beach. This term was initially missing.

**Impact**: Eroded sediment was "lost" from the system.

**Fix**: Added `deltavramp_ero` to the beach volume update:
```python
vbeachtotnew = vbeachtot_old + (-qwnet - qslr + qnour - qbar) * p.dt + vnourb[i] + deltavramp_ero
```

---

### 4. Initial Beach Volume Inconsistency
**Location**: `src/cs_model/main.py`, line 268

**Problem**: The initial `vbeach` value was taken from the config file, but this didn't necessarily match the beach geometry defined by the initial shoreline position.

**Impact**: The model started with an inconsistent state where the beach volume didn't match its geometry.

**Fix**: Calculate initial `vbeach` from the shoreline position:
```python
vbeach = ((p.yg - p.ys) - p.bl) / p.al  # Calculate from geometry instead of config
```

---

### 5. Numpy 2.x Compatibility
**Location**: `src/cs_model/functions.py`, line 87

**Problem**: `np.trapz` is deprecated in numpy 2.x and raises an AttributeError.

**Fix**: Added compatibility layer:
```python
try:
    return np.trapezoid(y_clip, x)
except AttributeError:
    return np.trapz(y_clip, x)
```

---

## Mass Conservation Results

### Test Simulation (24 hours, 144 timesteps)

**Before fixes:**
- Total mass balance error: 12.23 m³/m (2.48% of initial volume)

**After fixes:**
- Dune component: 0.0052% error
- Beach component: 0.0006% error
- Bar component: 0.4772% error
- **Total system: 0.0022% error** ✓

The small remaining error is due to numerical integration tolerance and is well within acceptable limits for a morphodynamic model.

---

## Where Was Sediment Being Lost/Created?

The main issue was in the **beach volume calculation**. Specifically:

1. **Missing qbar flux** meant that sediment deposited in the bar was not being subtracted from the beach, artificially inflating beach volume
2. **Geometric discontinuity** from using the updated ramp_toe position created artificial volume changes when the dune toe moved
3. **Combined effect**: The system was creating about 0.076 m³/m of sediment per timestep (roughly 11 m³/m over 144 steps)

---

## How to Verify Mass Conservation

### Run the Mass Conservation Report

```bash
python mass_conservation_report.py test_mass_conservation.yaml
```

This will:
1. Run the model with your configuration
2. Track all sediment volumes (dune, beach, bar) over time
3. Calculate integrated fluxes (aeolian, bar-beach, erosion)
4. Generate a comprehensive report showing mass balance errors
5. Verify that total system mass is conserved (< 0.1% threshold)

### Example Output

```
================================================================================
CS MODEL MASS CONSERVATION REPORT
================================================================================

INITIAL STATE:
  Dune volume:      333.3440 m³/m
  Beach volume:    4235.9769 m³/m
  Bar volume:        22.6079 m³/m
  Total volume:    4591.9288 m³/m

FINAL STATE:
  Dune volume:      333.7048 m³/m
  Beach volume:    4220.2435 m³/m
  Bar volume:        37.8793 m³/m
  Total volume:    4591.8276 m³/m

VOLUME CHANGES:
  Δ Total:      -0.1012 m³/m (  -0.00%)

MASS BALANCE ASSESSMENT:
  TOTAL ERROR:    -0.101236 m³/m (0.002205%)

✓ RESULT: Mass is conserved within acceptable limits!
```

---

## Files Modified

1. **src/cs_model/main.py**
   - Beach volume calculation timing (lines 470-503)
   - Beach volume update equation (line 507)
   - Initial beach volume calculation (line 268)

2. **src/cs_model/functions.py**
   - Numpy 2.x compatibility (line 87)

3. **New Files Added**
   - `mass_conservation_report.py` - Comprehensive mass conservation checking tool
   - `.gitignore` - Exclude cache files and output files
   - `test_mass_check.py`, `check_mass.py` - Test scripts for debugging

---

## Recommendations

1. **Always use the mass conservation report** when running simulations to verify results
2. **Monitor bar component**: The bar has the largest relative error (0.48%), though still small in absolute terms
3. **Long simulations**: For multi-year simulations, the 0.0022% error may accumulate. Consider periodic checks.
4. **Config file**: The `vbeach` parameter in the config file is now ignored (calculated from geometry). You can remove it or update documentation to note this.

---

## Questions?

If you have any questions about the mass conservation fixes or need help running the verification script, please let me know!
