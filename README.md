# CS-hybrid

This repository contains an extended **1D version of the CS Model** (Hallin et al., 2019), adapted to represent **vegetation dynamics** and **hard coastal structures** for **hybrid dune–dike systems**.

The model simulates cross-shore profile evolution under:

- Waves and water levels
- Sediment transport and morphological change
- Vegetation growth and feedbacks
- Hard structures (e.g., dikes and revetments)

---

## Authors

- **Filipe Galiforni Silva** (TU Delft)
- **Ad Reniers** (TU Delft)
- **Sierd de Vries** (TU Delft)
- **Floris van Rees** (Deltares)
- **Caroline Hallin** (Lund University)

TU Delft, 2026

---

## Repository Structure

```text
csmodel-hybrid/
├── src/cs_model/
│   ├── main.py                 # Main model driver
│   ├── functions.py            # Morphodynamics and utility routines
│   └── mass_conservation.py    # Sediment conservation diagnostics
├── examples/
│   ├── katwijk/
│   ├── oosteroever/
│   └── sediment_conservation/
├── analyze_conservation_timeseries.py
└── mass_conservation_report.py
```

---

## Environment Setup (Conda)

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate csmodel-hybrid
```

---

## Quick Start

Example runs:

```bash
cd examples/oosteroever
python run_example_oosteroever.py
```

or

```bash
cd examples/katwijk
python run_example_katwijk.py
```

Each run writes a pickle output file and can also generate sediment conservation diagnostics.

---

## Sediment Conservation (Updated)

Sediment conservation is tracked during the simulation and can be reported at the end of a run.

### 1) In-run conservation tracking

`run_model(...)` stores compartment-level conservation time series in:

- `results['conservation']['dune_error']`
- `results['conservation']['beach_error']`
- `results['conservation']['bar_error']`
- `results['conservation']['total_error']`
- and their percentage counterparts (`*_error_pct`)

The time series include the initial state (`t=0`) plus all simulated time steps.

### 2) End-of-run conservation report

When `check_mass_conservation=True`, the model runs
`cs_model.mass_conservation.check_mass_conservation(...)` and prints/saves a report that includes:

- Initial and final total sediment volume
- Integrated aeolian transport, erosion, bar exchange, and nourishment
- Phantom-beach correction linked to ramp-toe movement
- Drift diagnostics (per step, per hour, annualized)
- Optional timestep spike diagnostics if conservation fails

The full text report is saved next to the output pickle file as:

`<output_name>_mass_conservation.txt`

### 3) Current acceptance criteria

The conservation verdict uses a dual criterion:

- **Primary:** `|ΔV| / transport < 0.5%`
- **Secondary:** `|ΔV| / V0 < 1.0%`

A run is marked as conserved when **either** criterion is met.

### 4) Running with conservation enabled

```python
from cs_model.main import run_model

results = run_model(
    "config.yaml",
    check_mass_conservation=True,
    skip_first_step_conservation=False,
)
```

For focused regression checks of conservation behavior, see:

- `examples/sediment_conservation/`
- `examples/sediment_conservation/test16/run_example_test16.py`

---

## Reference

Hallin, C., Larson, M., & Hanson, H. (2019).  
*Simulating beach and dune evolution at decadal to centennial scale under rising sea levels.*  
PLoS ONE, 14(4), e0215651.  
https://doi.org/10.1371/journal.pone.0215651




