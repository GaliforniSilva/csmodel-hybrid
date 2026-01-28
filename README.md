# CS-hybrid

This repository contains an extended **1D version of the CS Model** originally developed by **Hallin et al. (2019)**, adapted to include **vegetation dynamics** and **hard coastal structures** in order to simulate **hybrid dune–dike systems**.

The model simulates cross-shore coastal profile evolution under the influence of:

- Waves and water levels  
- Sediment transport and morphological change  
- Vegetation growth and feedbacks  
- Hard structures (e.g., dikes, revetments)  

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
│
├── src/
│   ├── __init__.py
│   ├── main.py          # Main model driver
│   └── functions.py     # Helper and utility functions
│
└── examples/
    ├── oosteroever/
    └── katwijk/
```

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

## Prerequisites
- Python 3.11 or higher
- Required packages: numpy, scipy, matplotlib, shapely, pandas, netCDF4, pyyaml

## Quick Start

##### Running a Simulation with a YAML Config File

Each example uses a YAML configuration file specifying all necessary parameters and file paths. Paths in the YAML should be **relative to the YAML file location** for portability.

Example YAML files:

- `examples/katwijk/config_hardstruct.yaml`
- `examples/oosteroever/config_vegetation.yaml`

###### Run the simulation

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

After completion:
* The output pickle file will be saved in the directory specified in the YAML (e.g., `cs_model_oosteroever.pkl`).
* Each example automatically generates a plot for visualization.





