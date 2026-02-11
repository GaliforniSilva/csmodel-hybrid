# Bar Growth Investigation

## Problem Statement

User observed that approximately 755 m³/m of sediment goes to the bar over one year, which seems suspicious given that:
- Beach loses 777 m³/m total
- Dune gains only 24 m³/m via aeolian transport
- **97% of beach sediment goes to bar, only 3% to dune**

## Analysis

### Bar-Beach Exchange (qbar) Behavior

Investigation of `qbar` values over the year shows:
- **qbar is ALWAYS positive** (0 negative values out of 52,560 timesteps)
- Mean qbar: 0.000024 m³/m/s
- Max qbar: 0.00018 m³/m/s
- Integrated: 755.53 m³/m/year

This means sediment ALWAYS flows from beach to bar, never the reverse.

### Bar Volume Evolution

| Time      | Bar Volume (m³/m) | Change from Initial |
|-----------|-------------------|---------------------|
| Initial   | 22.6              | -                   |
| 1 month   | 374.0             | +351 (1555%)        |
| 6 months  | 760.7             | +738 (3267%)        |
| 1 year    | 778.0             | +755 (3341%)        |

The bar grows exponentially initially, then slows as it approaches equilibrium.

## Root Cause: CB Parameter

The bar equilibrium volume is calculated in `subaq()` function (functions.py:2107):

```python
VOLBAREQ = CB * DEAN**ACOEFF * STP**BCOEFF * LO**2
```

With current parameters:
- Wave height (RMS): 0.707 m
- Wave period: 6.0 s
- **CB = 2** (from config)
- Dean parameter: ~11.8
- **Equilibrium bar volume: 2129 m³/m**

Since `VOLBAREQ >> VOLBAR` throughout the simulation:
- `DVOL = (VOLBAREQ - VOLBAR) * (1 - exp(-LAMBDA * DT))` is always positive
- `QBAR = DVOL / DT` is always positive
- Sediment continuously flows from beach to bar

### CB Parameter Sensitivity

Effect of CB on equilibrium bar volume:

| CB Value | Equilibrium Bar Volume | Comment                          |
|----------|------------------------|----------------------------------|
| 0.05     | 53 m³/m                | Very small bar                   |
| **0.1**  | **106 m³/m**           | **Reasonable (recommended)**     |
| **0.2**  | **213 m³/m**           | **Reasonable (recommended)**     |
| 0.5      | 532 m³/m               | Moderate bar                     |
| 1.0      | 1064 m³/m              | Large bar                        |
| **2.0**  | **2129 m³/m**          | **Unrealistically large** ⚠️     |
| 5.0      | 5319 m³/m              | Extremely large                  |

## Why Wasn't This Noticed Before?

### Original Code Behavior

In the ORIGINAL code (main branch, before mass conservation fixes):

```python
# Line ~505 (original)
vbeachtotnew = vbeachtot + (-qwnet - qslr + qnour) * p.dt + vnourb[i] + deltavramp_ero
# NOTE: Missing -qbar term!
```

**The beach volume was NOT being reduced by qbar**, even though:
- The bar function was growing the bar: `volbar += qbar * dt`
- Sediment was being removed from the beach-bar system

This created a **mass balance error** where:
- Sediment "magically appeared" in the bar
- Beach volume was artificially inflated
- The excessive bar growth was HIDDEN

### After Mass Conservation Fix

My fix added the `-qbar * dt` term:

```python
# Line 509 (current)
vbeachtotnew = vbeachtot_old + (-qwnet - qslr + qnour - qbar) * p.dt + vnourb[i] + deltavramp_ero
```

This:
✅ **Fixed the mass conservation** (now 0.004% error)  
🔍 **Revealed the underlying CB parameter issue**

## Impact on Results

With CB=2, the excessive bar growth dominates the sediment budget:

| Flux                    | Volume (m³/m/year) | % of Beach Loss |
|-------------------------|--------------------|-----------------|
| Beach → Bar             | 755.5              | 97.1%           |
| Beach → Dune (aeolian)  | 24.4               | 3.1%            |
| Dune → Beach (erosion)  | -2.5               | (gain)          |
| **Net Beach Loss**      | **777.4**          | **100%**        |

This is physically unrealistic - most beach sediment should not be going offshore to the bar.

## Recommendations

### 1. Reduce CB Parameter ⭐

**Change in config file:**
```yaml
constants:
  cb: 0.1  # or 0.2 (instead of 2)
```

This will:
- Reduce equilibrium bar volume to ~100-200 m³/m
- Allow more balanced sediment distribution
- Maintain proper mass conservation

### 2. Verify Against Field Data

If available, check:
- Typical bar volumes for the study site
- Historical bar growth rates
- Equilibrium bar characteristics

### 3. Consider Seasonal Variations

The bar equilibrium might vary with:
- Storm events (larger waves → larger equilibrium)
- Calm periods (smaller waves → smaller equilibrium)

Consider implementing time-varying CB if needed.

## Testing with Corrected CB

To test with CB=0.1, update the config file:

```yaml
constants:
  # ... other parameters ...
  cb: 0.1  # Reduced from 2.0
  # ... other parameters ...
```

Expected results with CB=0.1:
- Equilibrium bar volume: ~106 m³/m
- Bar growth: ~23 → 80-100 m³/m (more reasonable)
- More sediment retained on beach
- Better balance between bar and aeolian transport

## Conclusion

The excessive bar growth (755 m³/m) is NOT an error in mass conservation, but rather:

1. ✅ **Mass conservation is correct** - my fix properly accounts for qbar
2. ⚠️ **CB parameter is too high** - CB=2 gives unrealistic bar equilibrium
3. 🔍 **Issue was hidden** - original code's mass balance error masked the problem
4. 💡 **Solution**: Reduce CB to 0.1-0.2 for realistic bar behavior

The `-qbar` term should be KEPT in the code for proper mass conservation.
