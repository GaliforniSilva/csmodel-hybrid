# Terminal Output Example

This document shows what users will see when running a simulation with `check_mass_conservation=True`.

## Complete Terminal Output Flow

```
Progress: 0.0% (0 of 144 time steps)
Progress: 0.7% (1000 of 144000 time steps)
...
Progress: 99.3% (143000 of 144000 time steps)
Simulation completed. Results saved to test_short.pkl.

================================================================================
CS MODEL MASS CONSERVATION REPORT
================================================================================

Simulation details:
  Number of time steps: 144
  Time step: 600 seconds
  Total simulation time: 24.0 hours (1.0 days)

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
  Δ Dune:       +0.3608 m³/m (  +0.11%)
  Δ Beach:     -15.7335 m³/m (  -0.37%)
  Δ Bar:       +15.2714 m³/m ( +67.55%)
  Δ Total:      -0.1012 m³/m (  -0.00%)

INTEGRATED FLUXES (over entire simulation):
  Aeolian transport from beach to dune:       0.3782 m³/m
    - Deposited on ramp:        0.0000 m³/m
    - Deposited on crest:       0.0936 m³/m
    - Deposited on lee:         0.2847 m³/m
  Bar-beach exchange (to bar if +):          15.3793 m³/m
  Wave erosion (from dune):                  -0.0000 m³/m

MASS BALANCE ASSESSMENT:

Expected changes (based on fluxes):
  Expected Δ Dune:       +0.3782 m³/m
  Expected Δ Beach:     -15.7575 m³/m
  Expected Δ Bar:       +15.3793 m³/m
  Expected Δ Total:      +0.0000 m³/m

Mass balance errors:
  Dune error:     -0.017402 m³/m (0.0052%)
  Beach error:    +0.024060 m³/m (0.0006%)
  Bar error:      -0.107894 m³/m (0.4772%)
  TOTAL ERROR:    -0.101236 m³/m (0.002205%)

✓ RESULT: Mass is conserved within acceptable limits!
  (Error < 0.1% threshold)

================================================================================
Mass conservation summary saved to test_short_mass_conservation.txt

================================================================================
MASS CONSERVATION SUMMARY
================================================================================
Status: ✓ PASSED
Total Volume Error: -0.101236 m³/m (0.002205%)
Initial Total Volume: 4591.9288 m³/m
Final Total Volume: 4591.8276 m³/m
Error Threshold: 0.1%

✓ Mass is conserved within acceptable limits!
================================================================================

```

## Key Features

1. **Detailed Report**: Comprehensive analysis with all volume changes, fluxes, and errors
2. **File Save Notification**: Confirms the summary file was saved
3. **Final Summary**: Quick reference with pass/fail status and key metrics

## Benefits

- **Visibility**: Users can immediately see if mass is conserved
- **Convenience**: Both detailed and summary information are available
- **Persistence**: Summary file is saved for record-keeping
- **Clarity**: Clear pass/fail status with visual indicators (✓/⚠)
