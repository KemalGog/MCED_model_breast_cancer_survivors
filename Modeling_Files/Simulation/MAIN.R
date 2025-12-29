# =====================================================================
# MAIN.R
#
# Author: Kemal Caglar Gogebakan
#
# Main file to project reductions in late-stage cancer incidence for breast cancer survivors 
# and age-matched average-risk women via multi-cancer screening
# =====================================================================

suppressPackageStartupMessages({
  library(here)
  library(fs)
})

options(stringsAsFactors = FALSE)
options(mced.verbose = TRUE)

# ---------------------------------------------------------------------
# 0) Source function files 
# ---------------------------------------------------------------------
source(here("Modeling_Files", "Simulation", "01. data_read.R"))            
source(here("Modeling_Files", "Simulation", "02. load_fitted_data.R"))
source(here("Modeling_Files", "Simulation", "03. mced_eligibility.R"))   
source(here("Modeling_Files", "Simulation", "04. stage_shift.R"))         
source(here("Modeling_Files", "Simulation", "05. single_cancer_screening.R"))
source(here("Modeling_Files", "Simulation", "06. multi_cancer_screening.R"))
source(here("Modeling_Files", "Simulation", "07. run_all_scenarios.R"))

# ---------------------------------------------------------------------
# 1) User-defined inputs
# ---------------------------------------------------------------------

# Breast cancer patients age group at diagnosis
bc_dx_start   <- 45
bc_dx_end     <- 59
years_to_mced <- 5      # followed for Y years, then enter MCED screening

# MCED screening program design
end_age    <- 74        # end of MCED screening
screen_int <- 1         # yearly screening

# MCED test sensitivities (assumed same across cancers)
early_sens_vec <- c(0.30, 0.50)   # early-stage sensitivities
late_sens      <- 0.95            # late-stage sensitivity 

# Natural history parameter pairs (OMST = overall mean sojourn time; DMST = late-stage mean sojourn time)
OMST_DMST_pairs <- list(
  list(OMST = 1, DMST = 0.5),
  list(OMST = 2, DMST = 1.0),
  list(OMST = 2, DMST = 0.5)
)

# Configuration for file paths, BC survivor risk multipliers, and list of cancer sites modeled for MCED screening
cfg <- get_config(
  SIR_0_5_pos = 1.12,    # standardized incidence ratio to inflate incidence of all cancers combined for survivors (HR+ survivors, used in 0–5 years)
  SIR_0_5_neg = 1.38,    # standardized incidence ratio to inflate incidence of all cancers combined for survivors (HR− survivors, used in 0–5 years)
  HR_pos_acm  = 1.79,    # all-cause mortality hazard ratio for survivors (HR+ survivors)
  HR_neg_acm  = 1.79,    # all-cause mortality hazard ratio for survivors (HR− survivors)
  data_dir    = here("Data"),
  fitted_dir  = here("Modeling_Files", "Fitting", "Fitted_SEER_Data"),
  cancer_list = c(
    "anus","bladder","cervix","colorectal","esophagus","head","kidney",
    "liver","lung","lymphoma","melanoma","ovary","pancreas","stomach",
    "thyroid","uterus"
  )
)

# ---------------------------------------------------------------------
# 2) Load and adjust data and generate inputs for the model (all-cause mortality, incidence of all-cancers, non-met BC incidence)
# ---------------------------------------------------------------------
inputs <- load_and_adjust_data(include_baseline = TRUE, cfg = cfg)

# Survivor all-cause mortality by HR status (HR+ and HR−)
qx_survivor_by_hr <- inputs$risk_adjusted$qx_mortality_by_hr

# Baseline all-cause mortality (average-risk population)
# Average-risk women do not have HR status, so we use the same table for HR+ and HR−.
qx_baseline <- inputs$baseline_risk$qx_mortality_by_hr[["HR+"]]

# List of cancer types included in MCED model
cancer_list <- inputs$cancer_list

# ---------------------------------------------------------------------
# 3) Compute MCED screening entry-age distribution for BC survivors
# ---------------------------------------------------------------------
# Uses BC non-metastatic incidence to construct age distribution at diagnosis by hormone-receptor status.
# Follows survivors for Y years and removes those who die of any cause or develop
# a new cancer.

elig <- build_mced_entry_ages(
  risk_adjusted         = inputs$risk_adjusted,
  diag_age_min          = bc_dx_start,
  diag_age_max          = bc_dx_end,
  years_to_screen_start = years_to_mced
)

# Entry-age distributions for HR+ and HR− survivors (will be used for age-matched average-risk women as well.)
mced_start_age_dist_by_hr <- elig$survivors_entry_age

# ---------------------------------------------------------------------
# 4) Load age- and stage-specific incidence data for all cancers to be screened in MCED screening program.
# ---------------------------------------------------------------------
fitted_data <- load_fitted_data(
  rdata_folder = inputs$rdata_folder,
  cancer_list  = cancer_list,
  verbose      = TRUE
)

# ---------------------------------------------------------------------
# 5) Define MCED screening start age group
# ---------------------------------------------------------------------
start_age_min <- bc_dx_start + years_to_mced
start_age_max <- bc_dx_end   + years_to_mced

cat("MCED entry ages: ", start_age_min, "to", start_age_max, "\n")
cat("Screening stops at age:", end_age, "\n\n")

# ---------------------------------------------------------------------
# 6) Run all scenarios: OMST/DMST × early sensitivity settings
# ---------------------------------------------------------------------
out_files <- run_all_scenarios(
  bc_dx_start                = bc_dx_start,
  bc_dx_end                  = bc_dx_end,
  years_to_mced              = years_to_mced,
  end_age                    = end_age,
  screen_int                 = screen_int,
  early_sens_vec             = early_sens_vec,
  late_sens                  = late_sens,
  OMST_DMST_pairs            = OMST_DMST_pairs,
  fitted_data                = fitted_data,
  qx_survivor_by_hr          = qx_survivor_by_hr,
  qx_baseline                = qx_baseline,
  mced_start_age_dist_by_hr  = mced_start_age_dist_by_hr,
  cancer_list                = cancer_list,
  date_obj                   = Sys.Date()
)

cat("Saved scenario files:\n")
cat(paste0(" • ", out_files), sep = "\n")
cat("\nDone.\n")
