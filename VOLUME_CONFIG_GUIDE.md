# Dune Volume Configuration Guide

## Overview

This guide explains how to correctly set the `vramp` parameter in your YAML configuration file. An incorrect `vramp` value will lead to inconsistent initial conditions and inaccurate sediment budget calculations.

## The Problem

The `vramp` parameter in `dune_volumes` **must** match the geometric parameters defined in `dune_parameters`. If these are inconsistent, the model will start with an incorrect sediment budget, leading to:

1. **Incorrect mass conservation**: The initial volume doesn't match the geometry
2. **Unrealistic dune evolution**: The model compensates for the inconsistency, creating artifacts
3. **Misleading results**: Growth/erosion rates will be distorted

## How to Calculate vramp

The ramp volume is calculated using a **simplified triangular approximation**:

```
vramp = 0.5 × (ys - ysprim)² × betas
```

Where:
- `ys` = X-position of ramp toe (seaward edge) [meters]
- `ysprim` = X-position of ramp crest (landward edge) [meters]
- `betas` = Ramp slope (dimensionless, dz/dx)
- `vramp` = Ramp volume per unit width [m³/m]

### Physical Interpretation

The formula represents the volume of a triangular cross-section:
- **Base length**: `(ys - ysprim)` meters
- **Height at toe**: `0` (at dfoot level)
- **Height at crest**: `betas × (ys - ysprim)` meters
- **Cross-sectional area**: `0.5 × base × height = 0.5 × (ys-ysprim)² × betas`

## Example Calculation

Given the following dune parameters:
```yaml
dune_parameters:
  ysprim: 130    # Ramp crest at x=130m
  ys: 190        # Ramp toe at x=190m
  betas: 0.1     # Ramp slope of 0.1 (1:10)
```

Calculate vramp:
```
Base length = ys - ysprim = 190 - 130 = 60 m
vramp = 0.5 × 60² × 0.1
      = 0.5 × 3600 × 0.1
      = 180.0 m³/m
```

Therefore, your configuration should be:
```yaml
dune_volumes:
  vramp: 180.0
```

## Common Mistakes

### ❌ Mistake 1: Using an arbitrary value
```yaml
dune_volumes:
  vramp: 80.0  # Wrong! This should be 180.0 for the geometry above
```

**Impact**: 55.6% error in initial volume → incorrect sediment budget

### ❌ Mistake 2: Not updating vramp when changing geometry
If you change any of these parameters:
- `ysprim`
- `ys`
- `betas`

You **must** recalculate and update `vramp`!

### ❌ Mistake 3: Confusing units
All parameters should be in SI units:
- Positions (`ys`, `ysprim`): meters [m]
- Slope (`betas`): dimensionless [m/m]
- Volume (`vramp`): cubic meters per meter width [m³/m]

## Validation Tool

Use the provided validation script to check your configuration:

```bash
python validate_yaml_volumes.py your_config.yaml
```

This will:
1. Calculate the expected `vramp` from your geometry
2. Compare it with your configured value
3. Report any discrepancies
4. Provide recommendations if incorrect

### Example Output (Incorrect Configuration)
```
✗ VALIDATION FAILED: Volume is inconsistent with geometry

RECOMMENDATION:
  Update dune_volumes.vramp from 80.00 to 180.00
```

### Example Output (Correct Configuration)
```
✓ VALIDATION PASSED: Volume is consistent with geometry
```

## Quick Reference Table

| Base Length (ys-ysprim) | Slope (betas) | vramp (m³/m) |
|------------------------|---------------|--------------|
| 40 m                   | 0.1           | 80.0         |
| 50 m                   | 0.1           | 125.0        |
| 60 m                   | 0.1           | 180.0        |
| 60 m                   | 0.2           | 360.0        |
| 70 m                   | 0.1           | 245.0        |
| 80 m                   | 0.1           | 320.0        |

## Integration with Model Code

The model uses `vramp` in two ways:

1. **Initial Condition** (main.py, line 279):
   ```python
   vramp = p.vramp  # Uses value from YAML
   ```

2. **Runtime Recalculation** (main.py, line 508):
   ```python
   vramp = partial_volume(dune_x, dune_z, substrate_x, substrate_z, 
                          dune_x['ramp_crest'], dune_x['ramp_toe'])
   ```

If the initial `p.vramp` is wrong, the model starts with inconsistent conditions, even though it recalculates `vramp` at each timestep.

## Best Practices

1. **Always validate** your YAML file before running simulations
2. **Document your calculation** in comments within the YAML file
3. **Recalculate vramp** whenever you modify geometry parameters
4. **Check existing configurations** when starting a new case
5. **Use the validation script** as part of your workflow

## Additional Notes

### Other Volume Parameters

While this guide focuses on `vramp`, similar geometric consistency principles apply to:
- `vbeach`: Should be consistent with shoreline position and beach profile
- `volbar`: Should match the subaqueous bar dimensions

However, these calculations are more complex and may require different validation approaches.

### When vramp Doesn't Match Geometry

If you intentionally want a different initial volume (e.g., for calibration), understand that:
1. The model will adjust positions/slopes in the first timestep to match
2. This creates an artificial transient at t=0
3. Mass conservation may show a "jump" at initialization
4. Document this choice clearly in your research notes

## References

For more details on the geometric calculations and mass conservation, see:
- `INCONSISTENCIES_QUICK_REFERENCE.md` - Known issues with volume calculations
- `DUNE_GROWTH_ANALYSIS.md` - Analysis of dune evolution
- `CODE_FIX_RECOMMENDATIONS.md` - Detailed technical recommendations

## Support

If you encounter issues or have questions:
1. Run the validation script first
2. Check the error messages carefully
3. Review this guide
4. Examine the example configurations in `examples/`
5. Report persistent issues with full details of your configuration
