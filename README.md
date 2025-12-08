# How to Run the Model

This project is organized into two stages:

1.  **Calibration of the natural history model** (required before simulation)
2.  **Simulation of multi-cancer early detection (MCED) screening**

Because fitted SEER outputs are **not** included in this repository, users should run the **calibration step** before running any MCED simulations.

------------------------------------------------------------------------

## 1. Download and Place Required Data Into the `Data/` Folder

This repository does **not** include SEER or Life Table data.

Users should download the following from official public sources:

-   SEER Incidence 2011-2015 (AJCC 6th Edition)
-   U.S. Life Tables (2018)

See `Data/README.md` for descriptions of required files.

Place all files into:

```         
Data/
```

using **the same filenames** described in the README.

------------------------------------------------------------------------

## 2. Run the Natural History Model Calibration

From R or RStudio:

``` r
source("Modeling_Files/Fitting/Fitting.R")
```

This script will:

-   read the SEER datasets from `Data/`
-   apply SIR multipliers for breast cancer survivors
-   fit the natural history model for 16 cancers × (HR+, HR--, average-risk)
-   save fitted objects into:

```         
Modeling_Files/Fitting/Fitted_SEER_Data/
```

------------------------------------------------------------------------

## 3. Run simulation model

Once model fitting is complete:

``` r
source("Modeling_Files/Simulation/MAIN.R")
```

The master script:

1.  loads fitted natural history model outputs
2.  generates MCED screening entry age distribution
3.  runs multi-cancer screening program
4.  produces all scenario results (OMST/LMST × sensitivity)

Outputs are stored under:

```         
Modeling_Files/Simulation/Batch_YYYYMMDD/
```

------------------------------------------------------------------------

## Questions / Contact

If you have issues running the model or questions about the methodology, please contact:

**Kemal Caglar Gogebakan**\
Email: [**kgogebak\@fredhutch.org**](mailto:kgogebak@fredhutch.org){.email}

------------------------------------------------------------------------
