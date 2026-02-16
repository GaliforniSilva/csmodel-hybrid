# Visual Guide: Main Geometric Artifacts in CS Model

## Issue #1: Beach Volume "Fake Loss" from Ramp Toe Movement

```
═══════════════════════════════════════════════════════════════════════
INITIAL STATE (t = 0)
═══════════════════════════════════════════════════════════════════════

                    Dune                          Beach
          |------------------------|----------------------------------|
          ^                        ^                                  ^
       lee_foot              ramp_toe=315m                     shoreline=415m
                                   
Beach Volume Formula:
V_beach_tot = ramp_toe × (d_clos + d_foot) + (shoreline - ramp_toe) × (d_clos + d_beach)
            = 315 × (8 + 3) + (415 - 315) × (8 + d_beach)
            = 3465 + 100 × (8 + d_beach)

═══════════════════════════════════════════════════════════════════════
AFTER RAMP GROWTH (t = dt)
Sediment deposited on ramp → ramp grows seaward → ramp_toe moves
═══════════════════════════════════════════════════════════════════════

                    Dune                      Beach
          |---------------------------|-----------------------------|
          ^                           ^                             ^
       lee_foot                ramp_toe=320m                 shoreline=415m
                                      (moved 5m!)
                                   
Beach Volume Formula NOW:
V_beach_tot = 320 × (8 + 3) + (415 - 320) × (8 + d_beach)
            = 3520 + 95 × (8 + d_beach)
            
DIFFERENCE in first term: 3520 - 3465 = +55 m³/m
BUT this volume didn't come from anywhere!

═══════════════════════════════════════════════════════════════════════
THE PROBLEM:
═══════════════════════════════════════════════════════════════════════

• The formula uses ramp_toe as the BOUNDARY of integration
• When ramp_toe moves, the INTEGRATION DOMAIN changes
• This creates a GEOMETRIC volume transfer
• But NO ACTUAL SEDIMENT moved from dune to beach!

In this example:
→ Beach "gains" 55 m³/m in the first term
→ Beach "loses" volume in second term (width decreased 100→95m)
→ Net effect depends on d_beach value
→ This is a BOOKKEEPING ARTIFACT, not real sediment transport

IMPACT: Can be 10-100 m³/m per event, cumulative over simulation
```

---

## Issue #2: Ramp Volume "Geometric Inflation"

```
═══════════════════════════════════════════════════════════════════════
SCENARIO: 1 m³/m of sediment deposited on ramp
═══════════════════════════════════════════════════════════════════════

BEFORE:
                                    crest
                                     /|
                                   /  |
                                 /    | h₀ = 5.25 m
                               /      |
                             /        |
                           /__________|
                       ramp_toe    crest_x
                          
Ramp height: h₀ = 5.25 m
Ramp slope: β = 0.107 (10.7%)
Ramp toe position: x₀

─────────────────────────────────────────────────────────────────────

AFTER (if height increases uniformly):
                                       crest_new
                                        /|
                                      /  |
                                    /    | h₁ = 6.25 m
                                  /      | (increased 1m)
                                /        |
                              /          |
                            /____________|
                      ramp_toe_new    crest_x
                      
New ramp height: h₁ = 6.25 m
Height increase: Δh = 1.0 m
Toe moves seaward by: Δx = Δh / β = 1.0 / 0.107 ≈ 9.35 m

═══════════════════════════════════════════════════════════════════════
GEOMETRIC VOLUME INCREASE:
═══════════════════════════════════════════════════════════════════════

The new ramp extends seaward, creating a triangular wedge:

                     |\
                     | \
                Δh   |  \
                     |___\
                      Δx

Volume of wedge = 0.5 × Δx × Δh
                = 0.5 × (Δh/β) × Δh
                = 0.5 × (Δh)² / β
                = 0.5 × (1.0)² / 0.107
                ≈ 4.67 m³/m

═══════════════════════════════════════════════════════════════════════
THE PROBLEM:
═══════════════════════════════════════════════════════════════════════

You deposited:    1.0 m³/m    (actual sediment)
Calculation gives: 4.67 m³/m   (geometric volume)
                   ─────────
Fake volume:      3.67 m³/m   (WHERE DID THIS COME FROM?!)

This is a GEOMETRIC ARTIFACT from extending the profile seaward.
The volume is QUADRATIC in the height change: V_fake ∝ (Δh)²

For steep slopes (small β):
→ Small height changes create LARGE fake volumes
→ β = 0.05: V_fake = 0.5 × 1² / 0.05 = 10 m³/m (!!)
→ β = 0.20: V_fake = 0.5 × 1² / 0.20 = 2.5 m³/m

IMPACT: Ramp appears to grow much faster than sediment input rate
```

---

## Issue #3: Wrong Slope Calculation

```
═══════════════════════════════════════════════════════════════════════
SCENARIO: Calculate slope of dune ramp for wind shear acceleration
═══════════════════════════════════════════════════════════════════════

PROFILE:
                              crest (elevation z = 8.25 m)
                               /|
                             /  |
                           /    |
                         /      |
                       /        |
                     /__________|  toe (elevation z = 3.0 m)
                                   
            |<-------- 49 m ------->|
                horizontal distance

═══════════════════════════════════════════════════════════════════════
CURRENT CALCULATION (WRONG):
═══════════════════════════════════════════════════════════════════════

Code: slope = dune_z['ramp_crest'] / dx_ramp
           = 8.25 / 49
           = 0.168  (16.8% slope)

PROBLEM: Uses absolute elevation (8.25 m) instead of elevation DIFFERENCE

This is NOT the slope of the ramp!
This is: (elevation above datum) / (horizontal distance)

═══════════════════════════════════════════════════════════════════════
CORRECT CALCULATION:
═══════════════════════════════════════════════════════════════════════

slope = (z_crest - z_toe) / dx
      = (8.25 - 3.0) / 49
      = 5.25 / 49
      = 0.107  (10.7% slope)

ERROR: 16.8% vs 10.7% → 57% overestimate of slope!

═══════════════════════════════════════════════════════════════════════
WHY THIS MATTERS:
═══════════════════════════════════════════════════════════════════════

Wind shear acceleration depends on slope:
→ Steeper slopes → faster wind → more acceleration
→ Wrong slope → wrong wind shear distribution
→ Wrong deposition pattern on dune

The model calculates:
τ_ramp = α × (I + β × slope) + 1

Where slope is used directly in the acceleration term.

With wrong slope:
→ Overpredicts acceleration on ramp
→ More sediment deposits on ramp than should
→ Less reaches crest and lee
→ Wrong dune shape evolution

IMPACT: Dune grows with wrong morphology, may be too steep or wrong profile
```

---

## Visual Summary: Where Fake Volumes Come From

```
┌─────────────────────────────────────────────────────────────────────────┐
│                     REAL vs FAKE VOLUME CHANGES                         │
└─────────────────────────────────────────────────────────────────────────┘

REAL VOLUME CHANGE (Physical sediment transport):
───────────────────────────────────────────────────
Source → Transport → Deposition/Erosion → Destination

Example: Aeolian transport
Beach ──qw×dt──> Dune
  -10 m³/m        +10 m³/m    ✓ Conservation: -10 + 10 = 0


FAKE VOLUME CHANGE (Geometric artifacts):
───────────────────────────────────────────────────
Boundary moves → Integration domain changes → Apparent volume change

Example: Ramp toe shifts
Beach calculation before: V = 315 × 11 + 100 × db = 3465 + ...
Beach calculation after:  V = 320 × 11 + 95 × db  = 3520 + ...
                                                      ─────
                                          Difference: +55 m³/m
                                          
✗ NO SEDIMENT MOVED! Just bookkeeping changed.


═══════════════════════════════════════════════════════════════════════
THE FUNDAMENTAL ISSUE:
═══════════════════════════════════════════════════════════════════════

Your model uses a MOVING BOUNDARY (ramp_toe) to:
1. Define dune extent
2. Calculate beach extent  
3. Integrate volumes

When this boundary moves:
→ Integration domains change
→ Volume formulas change
→ Creates apparent volume changes
→ Even when no sediment moved!

This is like moving a measuring tape while trying to measure:
- Measuring tape at 0-100m → length = 100m
- Measuring tape moves to 5-100m → length = 95m
- Did the object shrink? NO! The reference point moved.
```

---

## Conceptual Solution Approaches

```
═══════════════════════════════════════════════════════════════════════
APPROACH 1: Explicit Geometric Transfer Terms
═══════════════════════════════════════════════════════════════════════

Track boundary movements separately:

V_beach_new = V_beach_old 
            + flux_terms (qw, qbar, etc.)
            + explicit_erosion_from_dune
            - geometric_transfer_from_boundary_shift

Where geometric_transfer = f(Δramp_toe, profile_shape)

✓ Separates physical and geometric changes
✓ Maintains mass conservation
✗ Requires careful accounting


═══════════════════════════════════════════════════════════════════════
APPROACH 2: Fixed Reference Frame
═══════════════════════════════════════════════════════════════════════

Use fixed x-coordinates for volume calculation:

V_beach = integral from x_fixed_toe to x_shoreline

Track sediment DENSITY at each point instead of total volume.
When ramp grows, update density field.

✓ No geometric artifacts
✓ Clear physical meaning
✗ More complex implementation
✗ Requires profile discretization


═══════════════════════════════════════════════════════════════════════
APPROACH 3: Hybrid (Recommended for Your Model)
═══════════════════════════════════════════════════════════════════════

Keep pivot point concept but add explicit tracking:

1. Calculate volume change from fluxes: ΔV_flux
2. Calculate volume change from geometry: ΔV_geom  
3. Update profile to accommodate ΔV_flux
4. Track ΔV_geom as a bookkeeping transfer
5. Ensure ΔV_flux + ΔV_geom = ΔV_total

✓ Keeps your elegant pivot point concept
✓ Makes artifacts explicit and trackable
✓ Easier to validate and debug
✓ Can be added incrementally


═══════════════════════════════════════════════════════════════════════
KEY PRINCIPLE:
═══════════════════════════════════════════════════════════════════════

Geometric changes (boundary movements) should be:
1. TRACKED explicitly
2. SEPARATED from physical sediment transport
3. CONSERVED (what leaves one compartment enters another)
4. VALIDATED (sum to zero in closed system)

Your model is close - it just needs explicit handling of these geometric
transfers instead of mixing them with physical transport.
```

---

## Testing Checklist

```
Use these simple tests to verify fixes:

□ Zero-Flux Test
  • Set all transport rates to zero (qw=0, qbar=0, etc.)
  • Allow only geometric changes (ramp can grow from nourishment)
  • Check: Total system volume should remain constant
  • FAIL if: Volume changes without any flux

□ Pure Deposition Test  
  • Fixed deposition rate: deposit 10 m³/m per timestep
  • Run for 10 timesteps
  • Check: Total deposited = 10 × 10 = 100 m³/m
  • Check: Profile volume increase = 100 m³/m
  • FAIL if: Profile shows >100 m³/m (geometric inflation)

□ Erosion-Transfer Test
  • Cause dune erosion of 50 m³/m
  • Check: Dune volume decreases by 50 m³/m
  • Check: Beach volume increases by 50 m³/m
  • Check: Total volume unchanged
  • FAIL if: Beach increase ≠ dune decrease

□ Boundary Movement Test
  • Record initial ramp_toe position
  • Run simulation with moderate transport
  • Log every ramp_toe position change
  • Calculate geometric volume transfer from each movement
  • Check: Geometric transfers + flux transfers = total volume change
  • FAIL if: Residual > 0.1%

□ Slope Consistency Test
  • Calculate slope as: dz/dx (correct)
  • Calculate slope as: z/dx (current wrong method)
  • Compare values
  • Check: Correct method gives smaller, physically reasonable slopes
  • FAIL if: Both methods give same answer (means profile is wrong)
```

---

**For complete analysis, see MODEL_ANALYSIS_REPORT.md**
