# Quick Start: Fixing Your Dune Volume Configuration

## The Issue
Your YAML configuration has an incorrect `vramp` value that doesn't match your dune geometry.

## Quick Fix

**Step 1:** Update your YAML file

Change this:
```yaml
dune_volumes:
  vramp: 80.0
```

To this:
```yaml
dune_volumes:
  vramp: 180.0
```

**Step 2:** Verify the fix

```bash
python validate_yaml_volumes.py your_config.yaml
```

You should see:
```
✓ VALIDATION PASSED: Volume is consistent with geometry
```

## Why This Matters

The `vramp` value must match the formula:
```
vramp = 0.5 × (ys - ysprim)² × betas
     = 0.5 × (190 - 130)² × 0.1
     = 0.5 × 60² × 0.1
     = 180.0 m³/m
```

Your configuration had `vramp = 80.0`, which is only **44%** of the correct value. This causes:
- ❌ Incorrect initial sediment budget
- ❌ Artificial model behavior at t=0
- ❌ Mass conservation errors
- ❌ Wrong dune evolution results

## Complete Solution

Three new files are provided in this repository:

1. **`validate_yaml_volumes.py`** - Run this before any simulation to check your config
2. **`config_corrected.yaml`** - Your configuration with the correct vramp value
3. **`VOLUME_CONFIG_GUIDE.md`** - Full documentation on calculating vramp

## Next Steps

1. ✅ Update your YAML: `vramp: 80.0` → `vramp: 180.0`
2. ✅ Validate: `python validate_yaml_volumes.py your_config.yaml`
3. ✅ Re-run your simulations with the corrected configuration
4. ✅ Any previous results with `vramp=80.0` should be discarded

## Questions?

See the detailed documentation:
- **`VOLUME_CONFIG_GUIDE.md`** - How to calculate vramp correctly
- **`VOLUME_CONFIG_ISSUE_REPORT.md`** - Detailed analysis of the problem

---

**Important:** Always validate your YAML files before running simulations!
