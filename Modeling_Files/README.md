# Overview

This repository contains code supporting the analysis described in:

"Is Multi-Cancer Early Detection Likely to be More Beneficial in a Survivor Population? A Modeling Study in Breast Cancer Survivors" (Gogebakan et al.)

The project consists of two main components:

1.  Calibration of the natural history model
2.  Simulation of multi-cancer early detection (MCED) screening

------------------------------------------------------------------------

### 1. Calibration of natural history model

Located in:\
`Modeling Files/Fitting/`

**Files:**

`natural_history_code.R` --- Contains the core functions used to fit the natural history model to age- and stage-specific incidence data for a single cancer, given OMST and LMST inputs.

`Fitting.R` - Reads SEER incidence data for all 16 cancers included in the MCED model\
- Uses `natural_history_code.R` to fit the model for each cancer, for breast cancer survivors (HR+ / HR−) and average-risk women\
- Saves calibrated outputs into `Modeling Files/Fitting/Fitted SEER Data/`

------------------------------------------------------------------------

### 2. Simulation of multi-cancer early detection (MCED) screening

Located in:\
`Modeling Files/Simulation/`

**Files:**

`01.data_read.R`- Loads all required incidence and mortality inputs and adjusts them for breast cancer survivors using standardized incidence ratios (SIRs) and all-cause mortality hazard ratio.

`02.load_fitted_data.R`- Loads all fitted natural-history model outputs (metadata and rate objects) for survivors and average-risk women.

`03.mced_eligibility.R`- Projects start age distributions of MCED screening for survivors and age-matched average-risk women.

`04.stage_shift.R`-Implements the stage-shift model following Lange et al. (2024). *Reference:*\
Lange JM, Gogebakan KC, Gulati R, Etzioni R. *Cancer Epidemiology, Biomarkers & Prevention.* 2024.

`05.single-cancer_screening.R`-Runs the single-cancer screening model: for each cancer site, HR group, and start age, it projects late-stage incidence with and without MCED screening. Outputs are later aggregated in the multi-cancer module.

`06_multi-cancer_screening.R`- Aggregates single-cancer results across all modeled cancer sites using age-distribution weights.\
Also runs parallel projections for age-matched average-risk women.

`07_run_all_scenarios.R`- Runs all combinations of early sensitivity × OMST/DMST scenarios and saves the resulting `.rds` output files.

`MAIN.R`- Main script that executes the full modeling workflow end-to-end.
