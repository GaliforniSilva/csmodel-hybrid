# Complete Parameter Calibration Guide

## Executive Summary

Your suspicions were **100% correct**! The investigation revealed **TWO parameter issues**:

1. ✅ **CB = 2.0 is too high** → Causes 755 m³/m bar growth (unrealistic)
2. ✅ **Wind speed = 13.0 m/s is too high** → Causes 126 m³/m dune growth (unrealistic)

## Recommended Configuration

```yaml
constants:
  cb: 0.1  # Changed from 2.0

boundary_conditions:
  wind_speed: 10.0  # Changed from 13.0 m/s
```

## Results Comparison

### Original Configuration (CB=2, Wind=13)
```
❌ UNREALISTIC
- Bar: +755 m³/m (3341% increase!) 
- Dune: +22 m³/m
- Beach: -777 m³/m
- 97% of sediment goes to bar
```

### Interim Test (CB=0.1, Wind=13)
```
⚠️ PARTIALLY FIXED
- Bar: +16 m³/m ✓
- Dune: +126 m³/m (too high!)
- Beach: -145 m³/m
- 89% of sediment goes to dune (but too much!)
```

### Recommended Configuration (CB=0.1, Wind=10)
```
✅ REALISTIC
- Bar: +16 m³/m ✓
- Dune: +37 m³/m ✓ (matches your expected 40!)
- Beach: -54 m³/m ✓
- 70% to dune, 30% to bar (balanced)
```

## Visual Comparison

```
┌─────────────────────────────────────────────────────────────┐
│ DUNE GROWTH (m³/m/year)                                     │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│ CB=2, W=13:   ▓▓▓▓▓▓▓░░░░░░░░░░░ 22 m³/m                   │
│                                                             │
│ CB=0.1, W=13: ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 126 m³/m (!)     │
│                                                             │
│ CB=0.1, W=10: ▓▓▓▓▓▓▓▓▓░░░░░░░░░░░░░ 37 m³/m ✓             │
│                                                             │
│ Expected:     ▓▓▓▓▓▓▓▓▓▓░░░░░░░░░░░░ 40 m³/m               │
│                                                             │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│ BAR GROWTH (m³/m/year)                                      │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│ CB=2, W=13:   ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ 755 m³/m (!)  │
│                                                             │
│ CB=0.1, W=13: ▓░░░░░░░░░░░░░░░░░░░░ 16 m³/m ✓              │
│                                                             │
│ CB=0.1, W=10: ▓░░░░░░░░░░░░░░░░░░░░ 16 m³/m ✓              │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

## What Was Wrong?

### Problem #1: CB Parameter Too High

**Original: CB = 2.0**
- Creates bar equilibrium of 2129 m³/m (unrealistic)
- Bar continuously pulls sediment offshore
- Causes 755 m³/m bar growth per year
- Beach narrows significantly

**Fixed: CB = 0.1**
- Creates bar equilibrium of 106 m³/m (realistic)
- Moderate bar growth of 16 m³/m per year
- Beach remains stable

### Problem #2: Wind Speed Too High

**Original: Wind = 13.0 m/s**
- Used as constant for entire year (unrealistic)
- Equivalent to 47 km/h or 29 mph continuously
- Real coasts have calm periods and storms
- Overestimates annual transport by ~3x

**Fixed: Wind = 10.0 m/s**
- Effective "annual equivalent" wind speed
- Accounts for calm periods
- Gives realistic annual transport

## Why Were These Issues Hidden?

### Cascading Effects

The two issues interacted to mask each other:

1. **CB=2 consumed beach sediment**
   - 97% went to bar
   - Beach narrowed from 100m to 39m width

2. **Narrow beach limited aeolian transport**
   - Duran & Moore reduction: `qwnet *= (1 - 4.2*h/W)`
   - Narrow beach → large reduction
   - Only 24 m³/m transported despite high winds

3. **Mass balance error hid both**
   - Original code didn't subtract qbar from beach
   - Allowed unrealistic bar growth
   - Made results appear more normal

### After My Mass Conservation Fix

Fixing the mass balance revealed both issues:
1. Added `-qbar` term → revealed CB=2 problem
2. Fixed CB → revealed wind=13 problem

## Technical Details

### CB Parameter Effect

Bar equilibrium volume:
```
VOLBAREQ = CB * DEAN^1.333 * STP^1.0 * LO^2
```

| CB  | Equilibrium | Annual Growth | Realistic? |
|-----|-------------|---------------|------------|
| 0.1 | 106 m³/m    | ~16 m³/m      | ✅ Yes     |
| 0.2 | 213 m³/m    | ~30 m³/m      | ✅ Maybe   |
| 0.5 | 532 m³/m    | ~100 m³/m     | ⚠️ High    |
| 2.0 | 2129 m³/m   | ~755 m³/m     | ❌ No      |

### Wind Speed Effect

Aeolian transport approximately:
```
q_aeolian ∝ (wind_speed - threshold)³
```

| Wind (m/s) | Annual Aeolian | Realistic?       |
|------------|----------------|------------------|
| 8.5        | 12 m³/m        | Too low          |
| 10.0       | 37 m³/m        | ✅ **Good!**     |
| 10.5       | 49 m³/m        | Acceptable       |
| 13.0       | 128 m³/m       | Too high         |

Small wind changes → Large transport changes due to cubic scaling.

## Mass Conservation Verification

With recommended configuration (CB=0.1, Wind=10):

| Component | Error   | Status |
|-----------|---------|--------|
| Dune      | 0.001%  | ✅     |
| Beach     | 0.0000% | ✅     |
| Bar       | 0.01%   | ✅     |
| **Total** | **0.0004%** | ✅ **Excellent** |

## Sediment Budget Breakdown

### Original (CB=2, W=13) - WRONG
```
Beach loses 777 m³/m:
├─ 755 m³/m (97%) → Bar offshore    ❌
└─  22 m³/m (3%)  → Dune (aeolian)  ❌
```

### Recommended (CB=0.1, W=10) - CORRECT
```
Beach loses 54 m³/m:
├─ 16 m³/m (30%) → Bar offshore     ✅
└─ 37 m³/m (70%) → Dune (aeolian)   ✅
```

Much more balanced and realistic!

## Implementation Steps

### 1. Update Your Config File

**File: `your_config.yaml`**

```yaml
constants:
  # ... other parameters ...
  cb: 0.1  # CHANGED from 2.0
  # ... other parameters ...

boundary_conditions:
  # ... other parameters ...
  wind_speed: 10.0  # CHANGED from 13.0
  # ... other parameters ...
```

### 2. Re-run Your Simulations

```bash
python src/cs_model/main.py your_config.yaml
```

### 3. Verify Results

Use the mass conservation report tool:
```bash
python mass_conservation_report.py your_config.yaml
```

Expected results:
- Dune growth: 35-40 m³/m/year ✓
- Bar growth: 15-20 m³/m/year ✓
- Mass balance error: < 0.1% ✓

## Additional Recommendations

### Use Time-Varying Forcing (If Possible)

Instead of constants, use actual time series:
```yaml
boundary_conditions:
  wave_height: "path/to/wave_timeseries.nc"
  wind_speed: "path/to/wind_timeseries.nc"
```

This will:
- Account for calm vs stormy periods
- Give more realistic transport patterns
- Better represent seasonal variations

### Calibrate to Your Site

If you have field data:
1. **Measured dune growth rates** → Tune wind speed
2. **Bar volume surveys** → Tune CB parameter
3. **Beach profile evolution** → Validate overall model

### Sensitivity Testing

Test parameter ranges:
- CB: 0.05 to 0.2
- Wind: 9.0 to 11.0 m/s

To understand uncertainty and natural variability.

## Summary

✅ **Both parameters now calibrated**
- CB: 0.1 (fixes bar growth)
- Wind: 10.0 m/s (fixes dune growth)

✅ **Results match expectations**
- Dune: 37 m³/m/year (expected 40)
- Bar: 16 m³/m/year (reasonable)

✅ **Mass conservation perfect**
- Error: 0.0004%

✅ **Physically realistic**
- Balanced sediment distribution
- Sustainable beach evolution

## Questions Answered

**Q1: Was mass being conserved properly?**
- ✅ Yes, after fixes (0.0004% error)

**Q2: Why was bar growth so high (755 m³/m)?**
- ✅ CB=2.0 was too high

**Q3: Why is dune growth so high (126 m³/m) with CB=0.1?**
- ✅ Wind speed 13 m/s was too high

**Q4: What are the correct parameters?**
- ✅ CB=0.1, Wind=10 m/s → 37 m³/m dune growth

Your scientific intuition was spot-on in catching both issues! 🎯
