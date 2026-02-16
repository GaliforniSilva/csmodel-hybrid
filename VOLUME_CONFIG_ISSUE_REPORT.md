# Dune Volume Configuration Issue - Analysis Report

## Summary

This report documents the discovery and resolution of a critical configuration error in the dune volume parameters that was causing incorrect sediment budget calculations.

## Issue Identified

The YAML configuration file provided by the user contained a **geometrically inconsistent vramp value**:

```yaml
dune_parameters:
  ysprim: 130    # Ramp crest position (m)
  ys: 190        # Ramp toe position (m)
  betas: 0.1     # Ramp slope

dune_volumes:
  vramp: 80.0    # ❌ INCORRECT - should be 180.0
```

## Root Cause Analysis

### The Geometric Formula

The ramp volume must satisfy:
```
vramp = 0.5 × (ys - ysprim)² × betas
```

This formula represents a triangular profile:
- Base length: `ys - ysprim = 190 - 130 = 60 m`
- Height at crest: `betas × base = 0.1 × 60 = 6 m`
- Cross-sectional area: `0.5 × 60 × 6 = 180 m²`
- Volume per unit width: `180 m³/m`

### The Error

**Configured value**: 80.0 m³/m  
**Correct value**: 180.0 m³/m  
**Absolute error**: -100.0 m³/m  
**Relative error**: -55.6%

This means the model started with **only 44% of the sediment** that should be present in the ramp based on its geometry.

## Impact on Model Results

### 1. Initial Condition Mismatch

When the model initializes:
- The geometry (ysprim=130, ys=190, betas=0.1) defines a ramp profile
- The volume (vramp=80.0) doesn't match this geometry
- The model must reconcile this inconsistency, creating artificial transients

### 2. Mass Budget Errors

With incorrect initial volume:
- Total sediment budget starts from wrong baseline
- Beach-dune exchanges are calculated relative to this incorrect baseline
- Growth/erosion rates may be distorted to compensate

### 3. Sediment Conservation Issues

The model tracks volume changes via:
- **Geometric calculation**: Based on profile positions and slopes
- **Flux integration**: Based on transport rates

When initial vramp is wrong, these two approaches diverge at t=0, creating:
- Spurious volume "jumps" at initialization
- Confusion about whether mass is truly conserved
- Difficulty validating model physics

## Solution Implemented

### 1. Validation Script (`validate_yaml_volumes.py`)

Created a Python script that:
- Reads YAML configuration files
- Extracts geometric parameters (ysprim, ys, betas)
- Calculates expected vramp from geometry
- Compares with configured vramp value
- Reports discrepancies with clear recommendations

**Usage:**
```bash
python validate_yaml_volumes.py your_config.yaml
```

**Output for incorrect configuration:**
```
✗ VALIDATION FAILED: Volume is inconsistent with geometry

RECOMMENDATION:
  Update dune_volumes.vramp from 80.00 to 180.00
```

**Output for correct configuration:**
```
✓ VALIDATION PASSED: Volume is consistent with geometry
```

### 2. Corrected Configuration File

Created `config_corrected.yaml` with the fix:
```yaml
dune_volumes:
  vramp: 180.0  # ✓ CORRECT - matches geometry
```

### 3. Documentation

Created comprehensive documentation:
- **VOLUME_CONFIG_GUIDE.md**: User guide for understanding and calculating vramp
- **Comments in code**: Added explanations in main.py and functions.py
- **This report**: Detailed analysis of the issue and solution

### 4. Testing

Validated the solution:
- ✓ Original config (vramp=80.0) **FAILS** validation
- ✓ Corrected config (vramp=180.0) **PASSES** validation
- ✓ Existing example configs (Katwijk, Oosteroever) **PASS** validation

## Verification of Existing Configurations

Tested the validation script on existing example configurations:

### Katwijk Configuration
```
ysprim: 265.93 m
ys: 315.00 m
betas: 0.107
Expected vramp: 128.80 m³/m
Configured vramp: 128.80 m³/m
✓ VALIDATION PASSED
```

### Oosteroever Configuration
```
ysprim: 219.79 m
ys: 220.00 m
betas: 0.47
Expected vramp: 0.01 m³/m
Configured vramp: 0.01 m³/m
✓ VALIDATION PASSED
```

Both existing configurations are correct, confirming that this issue was specific to the user's configuration.

## Recommendations for Users

### Before Running Any Simulation:

1. **Always validate your YAML file:**
   ```bash
   python validate_yaml_volumes.py your_config.yaml
   ```

2. **If validation fails, update vramp:**
   - Use the recommended value from the validation report
   - Don't guess or use arbitrary values

3. **Document your geometry:**
   - Add comments in your YAML explaining how vramp was calculated
   - Include the formula: `vramp = 0.5 × (ys - ysprim)² × betas`

### When Changing Dune Geometry:

If you modify any of these parameters:
- `ysprim` (ramp crest position)
- `ys` (ramp toe position)
- `betas` (ramp slope)

You **must** recalculate and update `vramp` accordingly!

### Example YAML with Good Practices:

```yaml
dune_parameters:
  ysprim: 130
  ys: 190
  betas: 0.1

dune_volumes:
  # vramp calculated from geometry:
  # vramp = 0.5 × (190 - 130)² × 0.1 = 180.0
  vramp: 180.0
```

## Technical Details

### Why the Simple Formula Works

For a triangular ramp profile on flat substrate:
- The ramp is a straight line from (ysprim, sramp+dfoot) to (ys, dfoot)
- Height varies linearly: h(x) = betas × (ys - x)
- Volume integral: ∫[ysprim to ys] h(x) dx = 0.5 × (ys-ysprim)² × betas

The actual `ramp_volume()` function in the code performs numerical integration, which:
- Handles complex substrate topography
- Accounts for vertical offsets
- Provides higher accuracy

But for typical cases with flat substrate, the simple formula is very accurate.

### Relationship to Model Code

In `main.py` (line 279):
```python
vramp = p.vramp  # Initialize from config
```

Later (line 508):
```python
vramp = partial_volume(...)  # Recalculate from geometry
```

The model recalculates `vramp` at each timestep based on the actual geometry. However, starting with the wrong value creates:
- A discontinuity at t=0
- Potential mass balance artifacts
- Confusion in result interpretation

## Files Created/Modified

### New Files:
1. `validate_yaml_volumes.py` - Validation script
2. `config_corrected.yaml` - Corrected configuration
3. `config_original.yaml` - Original (incorrect) configuration for testing
4. `VOLUME_CONFIG_GUIDE.md` - User documentation
5. `VOLUME_CONFIG_ISSUE_REPORT.md` - This report

### Modified Files:
1. `src/cs_model/main.py` - Added explanatory comments
2. `src/cs_model/functions.py` - Enhanced docstring for ramp_volume()

## Conclusion

This analysis identified a **55.6% error in the initial ramp volume** that would have:
- Caused incorrect initial sediment budget
- Created artificial transients at model initialization
- Made mass conservation validation difficult
- Potentially affected all subsequent results

The solution provides:
- **Automated validation** to catch this error before running simulations
- **Clear documentation** explaining the correct calculation
- **Corrected configuration** ready to use
- **Code comments** to help future users understand the requirement

**Action Required**: Update your YAML configuration file to use `vramp: 180.0` instead of `vramp: 80.0`, then re-run your simulations with the corrected value.

## Questions for Follow-Up

1. **Origin of the 80.0 value**: Was this copied from another configuration, calculated with different parameters, or set arbitrarily?

2. **Existing results**: Any simulations run with vramp=80.0 should be considered invalid and re-run with the corrected value.

3. **Other parameters**: Should we also validate vbeach and volbar? These have more complex relationships but could benefit from similar checks.

4. **Automated testing**: Consider adding this validation as a pre-run check in the model itself, not just as a separate script.
