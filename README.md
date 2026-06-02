# CS-hybrid

This repository contains an extended **1D version of the CS Model** (Hallin et al., 2019), adapted to represent **vegetation dynamics** and **hard coastal structures** for **hybrid dune–dike systems**.

The model simulates cross-shore profile evolution under:

- Waves and water levels
- Sediment transport and morphological change
- Vegetation growth and feedbacks
- Hard structures (e.g., dikes and revetments)

This framework is intended for research applications at **annual to centennial time scales**, supporting studies on nature-based and hybrid coastal defence strategies.

---

## Authors

- **Filipe Galiforni Silva** (TU Delft)
- **Ad Reniers** (TU Delft)
- **Sierd de Vries** (TU Delft)
- **Floris van Rees** (Deltares)
- **Caroline Hallin** (Lund University)

TU Delft, 2026

---

## Model Description

The CS Hybrid Model is based on the original CS Model formulation presented in:

> Hallin, C., Larson, M., & Hanson, H. (2019).  
> *Simulating beach and dune evolution at decadal to centennial scale under rising sea levels.*  
> PLoS ONE, 14(4), e0215651.  
> https://doi.org/10.1371/journal.pone.0215651  

This version extends the model to explicitly represent:

- Vegetated dune evolution  
- Interaction between dunes and hard structures such as dikes  

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
```

---

## Main Scripts

1. `main.py`

The main executable file that runs the CS Hybrid simulation.

**Key features:**

- Reads model setup from a YAML configuration file  
- Simulates 1D cross-shore profile evolution  
- Includes vegetation and hard structure modules  
- Produces time-varying morphology and vegetation outputs  


2. `functions.py`

This file contains helper functions used throughout the model workflow, including:

- Geometric profile operations
- Sediment transport calculations
- Interpolation utilities
- Vegetation and structure-related routines

3. `mass_conservation.py`

This module provides mass conservation analysis for model simulations. It checks that sediment mass is conserved throughout the simulation by:

- Tracking volume changes in dune, beach, and bar compartments
- Comparing actual volume changes to expected changes based on sediment fluxes
- Generating detailed reports on mass balance errors
- Alerting users if mass balance errors exceed acceptable, user-defined thresholds

## Prerequisites
- Python 3.11 or higher
- Required packages: numpy, scipy, matplotlib, shapely, pandas, netCDF4, pyyaml
- Conda yaml available.


## Installation

The package can be installed directly with `pip`, which makes it importable from any Python environment without path manipulation:

```bash
# Editable install (recommended for development / coupling with other models)
pip install -e /path/to/csmodel-hybrid

# Or install from within the repository root
cd /path/to/csmodel-hybrid
pip install -e .
```

Once installed you can import the model anywhere:

```python
from cs_model import run_model

results = run_model("config.yaml")
```

Alternatively, install into a conda environment using the provided `environment.yml`:

```bash
conda env create -f environment.yml
conda activate csmodel-hybrid
pip install -e .
```



##### Running a Simulation with a YAML Config File

Each example uses a YAML configuration file specifying all necessary parameters and file paths. Paths in the YAML should be **relative to the YAML file location** for portability.

Example YAML files:

- `examples/katwijk/config_hardstruct.yaml`
- `examples/oosteroever/config_vegetation.yaml`

###### Run the simulation

After installing the package (see **Installation** above), import and run from anywhere:

```python
from cs_model import run_model

results = run_model("config_vegetation.yaml")
```

Or run an example script from its directory:

1. On a terminal, change directory to the desired example (katwijk or oosteroever):

```bash
cd PATH/examples/oosteroever
cd PATH/examples/katwijk

```

2. Then, run the script:

```bash
python run_example_katwijk.py
python run_example_oosteroever.py
```

### `run_model(...)` flags (runtime options)

`run_model` takes a YAML config path plus optional keyword flags:

- `check_mass_conservation=False`  
  Runs end-of-simulation sediment conservation diagnostics and writes a text report.
- `skip_first_step_conservation=False`  
  Excludes the first model step from conservation diagnostics when enabled.
- `plot_evolution=False`  
  Saves profile snapshots during the run as `<output_base>_<timestep>.png` inside
  `<output_base>_frames/`.
- `plot_evolution_interval=5.0`  
  Snapshot spacing as percent of total simulation duration.
- `save_evolution_movie=False`  
  Exports all saved snapshots to `<output_base>_evolution.gif` (requires `plot_evolution=True`).
- `save_evolution_video=False`  
  Exports all saved snapshots to `<output_base>_evolution.mp4` (requires `plot_evolution=True`).
- `manual_plot_limits=False`  
  Enables manual plotting limits for evolution snapshots.
- `manual_plot_y_lim_top=None`, `manual_plot_y_lim_bottom=None`  
  Manual vertical axis limits (m) when `manual_plot_limits=True`.
- `manual_plot_x_lim_right=None`, `manual_plot_x_lim_left=None`  
  Manual horizontal axis limits (m) when `manual_plot_limits=True`.

### Plotting and saving outputs

The output basename is taken from the YAML key:

```yaml
output:
  filename: cs_model_oosteroever.pkl
  output_disable: 0
  coupling: 0
```

- `output.filename`  
  Default pickle output path when saving is enabled.
- `output.output_disable`  
  When set to `1`/`true`, disables file outputs, including the pickle file and any derived saved outputs.
- `output.coupling`  
  When set to `1`/`true`, enables coupling mode: file outputs are disabled and general terminal printing is suppressed.

For this example, the model can produce:

- `cs_model_oosteroever.pkl`  
  Main serialized results dictionary.
- `cs_model_oosteroever_mass_conservation.txt`  
  Conservation summary report (when `check_mass_conservation=True`).
- `cs_model_oosteroever_<timestep>.png`  
  Evolution snapshots saved in `cs_model_oosteroever_frames/` (when `plot_evolution=True`).
- `cs_model_oosteroever_evolution.gif`  
  GIF movie from snapshots (when `save_evolution_movie=True`).
- `cs_model_oosteroever_evolution.mp4`  
  MP4 movie from snapshots (when `save_evolution_video=True`).

### Common run modes

1. Fast run without diagnostics or plots:

```python
run_model("config.yaml")
```

2. Conservation-focused run:

```python
run_model(
    "config.yaml",
    check_mass_conservation=True,
    skip_first_step_conservation=False,
)
```

3. Snapshot + GIF + MP4 export:

```python
run_model(
    "config.yaml",
    plot_evolution=True,
    plot_evolution_interval=5.0,
    save_evolution_movie=True,
    save_evolution_video=True,
)
```

After completion:
* The output pickle file is saved using `output.filename` from the YAML.
* Snapshot PNGs are saved in a dedicated `<output_base>_frames/` folder.
* Optional GIF/MP4 are saved next to the output pickle file.
* If enabled, the mass conservation report is printed and saved to a `_mass_conservation.txt` file.


## Sediment Conservation

The model includes automatic mass conservation checking that can be enabled when running simulations. This feature:

- **Validates** that sediment mass is conserved throughout the simulation
- **Tracks** volume changes in all compartments (dune, beach, bar)
- **Compares** actual changes to expected changes based on sediment fluxes
- **Reports** detailed mass balance information including:
  - Initial and final volumes for each compartment
  - Volume changes over the simulation
  - Integrated sediment fluxes (aeolian transport, wave erosion, bar-beach exchange)
  - Mass balance errors with threshold-based assessment
- **Saves** a summary report to a text file for record-keeping

### 1) In-run conservation tracking

`run_model(...)` stores compartment-level conservation time series in:

- `results['conservation']['dune_error']`
- `results['conservation']['beach_error']`
- `results['conservation']['bar_error']`
- `results['conservation']['total_error']`
- and their percentage counterparts (`*_error_pct`)

The time series includes the initial state (`t=0`) plus all simulated time steps.

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
from cs_model import run_model

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

