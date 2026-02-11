# Comparison: CB=2.0 vs CB=0.1

## Summary

Testing shows that reducing CB from 2.0 to 0.1 produces **much more realistic results**.

## Year Simulation Results Comparison

### Volume Changes

| Component | CB=2.0 (original) | CB=0.1 (corrected) | Change           |
|-----------|-------------------|--------------------| -----------------|
| Dune      | +22 m³/m (6.6%)   | +126 m³/m (37.8%)  | **+104 m³/m** ✅ |
| Beach     | -777 m³/m (-18%)  | -145 m³/m (-3.4%)  | **+632 m³/m** ✅ |
| Bar       | +755 m³/m (3341%) | +16 m³/m (73%)     | **-739 m³/m** ✅ |
| **Total** | -0.1 m³/m         | -2.4 m³/m          | -2.3 m³/m        |

### Sediment Flux Distribution

**With CB=2.0 (original):**
```
Beach loses 777 m³/m:
├─ 755 m³/m (97%) → Bar      ⚠️ Unrealistic!
└─  24 m³/m (3%)  → Dune
```

**With CB=0.1 (corrected):**
```
Beach loses 145 m³/m:
├─ 128 m³/m (89%) → Dune     ✅ Realistic!
└─  16 m³/m (11%) → Bar      ✅ Reasonable!
```

### Bar Equilibrium

| Parameter      | CB=2.0      | CB=0.1     |
|----------------|-------------|------------|
| Equilibrium    | 2129 m³/m   | 106 m³/m   |
| Initial        | 23 m³/m     | 23 m³/m    |
| Final          | 778 m³/m    | 39 m³/m    |
| Still growing? | **Yes** ⚠️  | Approaching |

### Aeolian Transport

| Metric              | CB=2.0       | CB=0.1        | Difference      |
|---------------------|--------------|---------------|-----------------|
| Total aeolian       | 24.4 m³/m    | 128.5 m³/m    | **+104 m³/m** ✅|
| Ramp deposition     | 0.0 m³/m     | 1.5 m³/m      | +1.5 m³/m       |
| Crest deposition    | 13.7 m³/m    | 67.3 m³/m     | +53.6 m³/m      |
| Lee deposition      | 10.7 m³/m    | 59.6 m³/m     | +48.9 m³/m      |

## Why Does CB Affect Aeolian Transport?

The beach volume affects aeolian transport through:

1. **Beach width changes** - With CB=2.0, the beach loses so much volume that it becomes narrower, reducing the fetch distance for aeolian transport

2. **Duran & Moore reduction** - The reduction factor depends on beach width:
   ```python
   qwnet = qwnet * max(1.0 - 4.2 * (dune_height) / (beach_width), 0.0)
   ```
   Narrower beach → larger reduction → less aeolian transport

3. **Feedback loop** - Less bar growth → wider beach → more aeolian transport → more dune growth

## Physical Interpretation

### CB=2.0 Behavior (Unrealistic)
- Bar grows aggressively, pulling sediment offshore
- Beach narrows significantly (400m → 339m shoreline)
- Reduced beach width limits aeolian transport
- Dune growth is suppressed
- System dominated by bar accretion

### CB=0.1 Behavior (Realistic)
- Bar grows moderately to near-equilibrium
- Beach width remains relatively stable
- Aeolian transport can operate effectively
- Dune grows substantially (~126 m³/m/year)
- Balanced sediment distribution

## Dune Growth Expectation

User expected approximately **40 m³/m/year** dune growth.

Results:
- With CB=2.0: **22 m³/m** (55% of expected) ❌
- With CB=0.1: **126 m³/m** (315% of expected) ⚠️

The CB=0.1 result is now HIGHER than expected! This suggests:

1. **CB=0.1 might be slightly too low** - try CB=0.15 or CB=0.2
2. **Duran & Moore reduction** - The ~17% reduction discussed earlier is still active
3. **Wave conditions** - Constant wave height of 1.0m may be unrealistic for annual average
4. **Wind conditions** - Constant wind speed of 13 m/s may be too high for annual average

## Recommendations

### 1. Fine-tune CB Parameter

Try intermediate values:
```yaml
cb: 0.15  # or 0.2
```

Expected results:
- CB=0.15 → Bar equilibrium ~160 m³/m
- CB=0.2 → Bar equilibrium ~213 m³/m

This might give dune growth closer to the expected 40 m³/m.

### 2. Validate Against Site Data

Check if available:
- Measured bar volumes
- Annual dune growth rates
- Beach width evolution
- Sediment budget studies

### 3. Consider Variable Forcing

Instead of constant conditions, use realistic time series:
- Seasonal wave climate
- Storm events
- Wind variability

This will likely reduce the annual-average aeolian transport.

## Conclusion

✅ **CB=0.1 produces physically realistic sediment distribution**
- Bar growth is moderate (~16 m³/m)
- Aeolian transport dominates (89% of beach loss)
- Dune growth is substantial

⚠️ **CB=0.1 may be slightly too aggressive for aeolian**
- Produces 126 m³/m dune growth (vs expected 40 m³/m)
- Consider CB=0.15-0.2 or using variable forcing conditions

🔍 **The original CB=2.0 was definitely too high**
- Caused unrealistic bar accretion (755 m³/m)
- Suppressed dune growth through beach narrowing
- Mass conservation fix revealed this issue

**Recommended next step**: Test CB=0.15 or CB=0.2 to find optimal balance.
