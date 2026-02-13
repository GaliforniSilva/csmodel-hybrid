# Mass Conservation Integration - Implementation Summary

## Overview

This document describes the implementation of automatic mass conservation checking in the CS Model simulation routine.

## Problem Statement

Previously, a mass conservation check routine was available but required manual execution after simulations. Users wanted this embedded into the running routine to automatically provide reports for all simulations done.

## Solution

### 1. New Reusable Module (`src/cs_model/mass_conservation.py`)

Created a new module with reusable functions for mass conservation analysis:

- **`check_mass_conservation(results, config, verbose=True)`**: Main function that analyzes simulation results and generates a detailed report
- **`calc_dune_vol()`**: Helper function to calculate dune volume
- **`calculate_actual_beach_volume()`**: Helper function to calculate beach volume

The module returns a dictionary with:
- `conserved`: Boolean indicating if mass is conserved within threshold
- `total_error_percent`: Total mass balance error as percentage
- `total_error`: Absolute total mass balance error
- Error metrics for individual compartments (dune, beach, bar)

### 2. Enhanced `run_model()` Function

Modified `src/cs_model/main.py` to add an optional parameter:

```python
def run_model(config_path, check_mass_conservation=False):
    """
    Run the CS model using a YAML configuration file.

    Parameters
    ----------
    config_path : str
        Path to the YAML configuration file
    check_mass_conservation : bool, optional
        If True, perform mass conservation check after simulation and print report.
        Default is False.

    Returns
    -------
    results : dict
        Model output
    """
```

When `check_mass_conservation=True`:
1. Simulation runs normally
2. After saving results, the mass conservation check is performed
3. A detailed report is printed to the console (verbose output)
4. A detailed report file is saved (e.g., `output_mass_conservation.txt`)
5. A final summary is printed to the terminal for quick reference

### 3. Updated Example Scripts

Both example scripts now enable mass conservation reporting by default:

**examples/oosteroever/run_example_oosteroever.py:**
```python
results = run_model("config_vegetation.yaml", check_mass_conservation=True)
```

**examples/katwijk/run_example_katwijk.py:**
```python
results = run_model("config_hardstruct.yaml", check_mass_conservation=True)
```

### 4. Refactored Existing Scripts

**mass_conservation_report.py** was refactored to use the new module, eliminating code duplication.

## Usage

### Basic Usage

```python
from cs_model.main import run_model

# Run simulation without mass conservation check (default behavior)
results = run_model('config.yaml')

# Run simulation with mass conservation check
results = run_model('config.yaml', check_mass_conservation=True)
```

### Advanced Usage

For programmatic access to mass conservation metrics:

```python
from cs_model.main import run_model
from cs_model.mass_conservation import check_mass_conservation
import yaml

# Run simulation
results = run_model('config.yaml', check_mass_conservation=False)

# Load config
with open('config.yaml', 'r') as f:
    config = yaml.safe_load(f)

# Perform mass conservation check with custom verbosity
mass_check = check_mass_conservation(results, config, verbose=True)

# Access metrics
if mass_check['conserved']:
    print(f"Mass is conserved! Error: {mass_check['total_error_percent']:.4f}%")
else:
    print(f"Warning! Error: {mass_check['total_error_percent']:.4f}%")
```

## Output

### Console Output

When mass conservation check is enabled, the following is printed to the terminal:

1. **Detailed Report**: A comprehensive analysis of mass conservation including simulation details, initial/final states, volume changes, integrated fluxes, and mass balance errors.

2. **Final Summary**: A concise summary printed at the end for quick reference, showing the pass/fail status and key metrics.

#### Detailed Report

```
================================================================================
CS MODEL MASS CONSERVATION REPORT
================================================================================

Simulation details:
  Number of time steps: 144
  Time step: 600 seconds
  Total simulation time: 23.8 hours (1.0 days)

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

Mass conservation report saved to test_short_mass_conservation.txt

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

### Report File

A detailed report file is automatically created (e.g., `output_mass_conservation.txt`) containing the full mass conservation analysis. The file contains the same detailed report shown in the terminal, including:

- Simulation details (time steps, duration)
- Initial and final volumes for all compartments (dune, beach, bar)
- Volume changes and percentages
- Integrated sediment fluxes (aeolian transport, bar-beach exchange, wave erosion)
- Mass balance assessment with expected changes
- Detailed error analysis for each compartment
- Final result (pass/fail) with threshold comparison

This makes it easy to review the complete mass conservation analysis later without re-running the simulation.

## Mass Conservation Thresholds

- **Threshold**: 0.1% of total initial volume
- **Conserved**: Total error < 0.1%
- **Warning**: Total error ≥ 0.1%

The 0.1% threshold is based on acceptable numerical error for long-term morphological simulations.

## Benefits

1. **Automatic Validation**: Every simulation can now be automatically validated for mass conservation
2. **Early Detection**: Issues are detected immediately after simulation completes
3. **Traceable**: Summary files provide permanent records of validation results
4. **Non-intrusive**: Existing code continues to work (parameter defaults to False)
5. **Flexible**: Can be used programmatically or via command-line scripts
6. **Consistent**: All simulations use the same validation methodology

## Testing

The implementation has been tested with:
- Short test simulations (24 hours)
- The test script `test_mass_integration.py` verifies the integration
- Both example configurations (Katwijk and Oosteroever) are ready to use

## Backward Compatibility

All changes are backward compatible:
- The `check_mass_conservation` parameter defaults to `False`
- Existing code that calls `run_model(config_path)` continues to work unchanged
- Users can opt-in to mass conservation reporting by adding the parameter

## Files Modified

1. `src/cs_model/main.py` - Enhanced run_model() function
2. `src/cs_model/mass_conservation.py` - New module (created)
3. `examples/katwijk/run_example_katwijk.py` - Enabled mass check
4. `examples/oosteroever/run_example_oosteroever.py` - Enabled mass check
5. `mass_conservation_report.py` - Refactored to use new module
6. `README.md` - Added documentation
7. `.gitignore` - Added output files
8. `test_mass_integration.py` - Test script (created)

## Future Enhancements

Potential future improvements could include:
- Configurable error thresholds
- Time-series plots of mass balance errors
- Export of detailed metrics to CSV/JSON
- Integration with continuous integration testing
