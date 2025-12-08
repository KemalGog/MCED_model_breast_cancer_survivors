## Fitted SEER Data (Structure Only)

This folder stores the **structure** of fitted natural history model outputs for all **16 cancers**, calibrated for:

-   HR-positive breast cancer survivors (`hr_pos`)
-   HR-negative breast cancer survivors (`hr_neg`)
-   Average-risk women (`hr_base`)

The actual `.Rdata` files are **not included**.\
This README describes the expected structure for modeling.

------------------------------------------------------------------------

### File Naming Convention

```         
{cancer}outfile_{group}.Rdata
```

Examples:

```         
lungoutfile_hr_pos.Rdata
ovaryoutfile_hr_neg.Rdata
kidneyoutfile_hr_base.Rdata
```

------------------------------------------------------------------------

### Contents of Each File

Each `.Rdata` file contains two objects:

#### 1. `<cancer>out_<group>`

-   A list of **three fitted models**, corresponding to OMST/LMST pairs: `(1, 0.5)`, `(2, 0.5)`, `(2, 1)`
-   Each model includes fitted transition rates and fitted age/stage-specific incidence.

#### 2. `metadata_<group>_<cancer>`

-   Contains OMST and LMST values for the corresponding fits.

------------------------------------------------------------------------

### Cancer Sites Included (16)

```         
anus, bladder, breast, cervix, colorectal, esophagus, head,
kidney, liver, lung, lymphoma, melanoma, ovary,
pancreas, stomach, thyroid, uterus
```

(Note: Breast is fitted but not used in MCED modeling.)

------------------------------------------------------------------------

### Use in Simulation

These files are automatically loaded by:

```         
Modeling Files/Simulation/02.load_fitted_data.R
```
