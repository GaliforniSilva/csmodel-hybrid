# Quick Reference: Model Inconsistencies by Priority

## HIGH PRIORITY - Fix These First

### 1. Beach Volume Fake Changes (§2.1)
**Location**: `main.py` lines 483-525  
**Problem**: When `ramp_toe` moves, beach volume calculation changes integration domain  
**Example**: Ramp grows 5m seaward → Beach "loses" 55 m³/m with no actual sediment movement  
**Impact**: Tens to hundreds of m³/m fake volume changes  

```python
# Problem area:
vbeachtot_old = dune_x['ramp_toe'] * (p.dclos + p.dfoot) + ...
# When ramp_toe changes, this term changes → fake volume transfer
```

### 2. Ramp Volume Geometric Inflation (§2.2)
**Location**: `functions.py` lines 1279-1505  
**Problem**: Ramp height increase creates volume proportional to (Δh)²/betas  
**Example**: 1m height increase → 4.7 m³/m volume for only 1 m³/m deposited sediment  
**Impact**: Square relationship creates large fake volumes  

```python
# Problem: geometric expansion when height increases
Fake Volume ≈ 0.5 × (Δh)² / betas
```

### 3. Wrong Slope in Deposition Calculation (§2.8)
**Location**: `functions.py` lines 876, 878  
**Problem**: Uses absolute elevation instead of elevation difference  
**Example**: `8.25/49 = 0.168` instead of `(8.25-3.0)/49 = 0.107`  
**Impact**: Wrong wind shear distribution → wrong deposition pattern  

```python
# Wrong:
slope_ramp = dune_z['ramp_crest'] / dx_ramp

# Should be:
slope_ramp = (dune_z['ramp_crest'] - dune_z['ramp_toe']) / dx_ramp
```

---

## MEDIUM-HIGH PRIORITY

### 4. Double Counting Risk (§2.4)
**Location**: `main.py` lines 484-486, 517, 521  
**Problem**: Eroded sediment added both geometrically (boundary shift) and explicitly (deltavramp_ero)  
**Impact**: Beach might receive sediment twice  

### 5. Erosion Edge Cases (§2.11)
**Location**: `functions.py` lines 1020-1260  
**Problem**: Complex conditional logic with >7 branches, potential gaps  
**Impact**: Failures in extreme scenarios  

### 6. Ramp-Dike Interaction (§2.12)
**Location**: `functions.py` lines 1291-1406  
**Problem**: Complex geometry assumptions may not match reality  
**Impact**: Incorrect volume near dike  

---

## MEDIUM PRIORITY

### 7. Fetch Distance Artifact (§2.7)
**Location**: `main.py` lines 448-450  
**Problem**: Fetch changes when ramp_toe moves, creating artificial feedback  
**Impact**: Aeolian transport affected by geometric artifact  

### 8. Dike Fill Volume (§2.3)
**Location**: `functions.py` lines 1294, 1325  
**Problem**: Formula assumes triangular difference, unclear physical meaning  
**Impact**: Over/underestimate fill requirements  

### 9. Foreshore Slope (§2.6)
**Location**: `main.py` line 350  
**Problem**: Calculates elevation/width instead of true slope  
**Impact**: Incorrect runup (works by coincidence in many cases)  

### 10. Shoreline Formula (§2.5)
**Location**: `main.py` line 521  
**Problem**: No validation for extreme cases  
**Impact**: Could create unphysical negative beach width  

---

## LOW PRIORITY (Documentation/Quality)

### 11. Wind Shear Logarithms (§2.9)
**Location**: `functions.py` lines 888-895  
**Impact**: Fragile, hard to verify  

### 12. Overwash Formula (§2.10)
**Location**: `main.py` lines 378-389  
**Impact**: May underestimate extreme overwash  

### 13. Volume Functions (§2.13)
**Location**: `functions.py` lines 98-118  
**Impact**: No input validation  

### 14. Intersection Comment (§2.14)
**Location**: `functions.py` lines 1232-1251  
**Impact**: Misleading comment  

---

## Testing Recommendations

### Immediate Tests:
1. **Zero-flux geometry test**: No transport, only boundary moves → Total volume should be constant
2. **Ramp growth test**: Known deposition → Check geometric vs actual volume ratio
3. **Erosion conservation test**: Dune loss should equal beach gain

### Validation Metrics:
- Per-timestep mass balance error
- Cumulative mass drift
- Geometric vs flux-based volume changes
- Position jump detection (flag if >10% of domain)

---

## Key Equations to Review

### Beach Volume:
```
vbeachtot = ramp_toe × (dclos + dfoot) + (shoreline - ramp_toe) × (dclos + db)
```
**Issue**: First term changes when ramp_toe moves

### Ramp Volume with Height Increase:
```
ΔV_geometric ≈ 0.5 × (Δh)² / betas
```
**Issue**: Quadratic in height change

### Slope Calculations:
```
# Current (wrong for some):
slope = z / dx

# Should be:
slope = (z_top - z_bottom) / dx
```

---

## Quick Visual: Fake Volume Example

```
Initial State:
Beach: |------- 100m -------| 
       ramp_toe          shoreline
Volume: 150 m³/m

After ramp growth (ramp_toe moves 5m seaward):
Beach: |------ 95m ------| 
       new_ramp_toe   shoreline
       
Calculation shows volume loss of ~55 m³/m
But no sediment actually moved!
This is a FAKE LOSS - just accounting boundary shift
```

---

## What to Monitor in Mass Conservation

1. **Total System Volume**
   - Should equal: Initial + Sources - Sinks
   - Sources: Nourishment, bar contribution
   - Sinks: Aeolian export, SLR offshore transport

2. **Compartment Transfers**
   - Dune → Beach (erosion)
   - Beach → Dune (aeolian)
   - Beach ↔ Bar (wave-driven)
   - **Check**: Transfers balance (what leaves one arrives at the other)

3. **Geometric Artifacts**
   - Track boundary position changes
   - Compare volume changes from:
     - Direct integration
     - Flux integration
     - Should match within numerical error

---

## Contact Points in Code

### Where volumes are calculated:
- `ramp_volume()` - functions.py:98
- `lee_volume()` - functions.py:119
- `total_volume()` - functions.py:142
- `partial_volume()` - functions.py:162

### Where boundaries move:
- `morpho_update()` - functions.py:997
- Beach volume calc - main.py:483-525

### Where mass enters/leaves:
- Aeolian: main.py:433-461
- Wave erosion: main.py:367-390
- Bar exchange: functions.py:2121-2160
- Nourishment: functions.py:2214-2229

---

**For full details, see: MODEL_ANALYSIS_REPORT.md**
