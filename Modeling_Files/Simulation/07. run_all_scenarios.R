##########################################################################
# 07. run_all_scenarios.R
#
# Author: Kemal Caglar Gogebakan
#
# Purpose:
#   Batch execution for MCED screening scenarios for breast cancer
#   survivors and age-matched average-risk women.
#
#   This script:
#       (1) Computes the MCED screening entry ages from diagnosis year window
#       (2) Runs MCED screening for:
#              – HR-positive survivors
#              – HR-positive average-risk (baseline)
#              – HR-negative survivors
#              – HR-negative average-risk (baseline)
#       (3) Saves each scenario as **one flat RDS file** in a single directory
#           (no nested scenario folders)
#       (4) Outputs exactly six files for 3 OMST–DMST pairs × 2 early sensitivities
#
# Inputs:
#   bc_dx_start, bc_dx_end      : Breast cancer diagnosis-year window
#   years_to_mced               : Years after diagnosis before MCED begins
#   end_age                     : Screening stops before this age
#   early_sens_vec              : Vector of early sensitivities (e.g., c(0.30, 0.50))
#   late_sens                   : Late sensitivity (scalar)
#   OMST_DMST_pairs             : List of natural history parameter pairs
#   fitted_data                 : Output from load_fitted_data()
#   qx_survivor_by_hr           : ACM tables for survivors ("HR+", "HR-")
#   qx_baseline                 : ACM for average-risk women
#   mced_start_age_dist_by_hr   : Entry-age distributions for HR+ and HR– survivors
#                                 (also used for age-matched baseline)
#   cancer_list                 : List of cancer sites evaluated
#
# Output:
#   • Creates directory: BC_Survivors_Batch_YYYYMMDD
#   • Produces 6 `.rds` files — one per scenario
#   • Each RDS stores:
#         meta info,
#         survivors HR+,
#         baseline HR+,
#         survivors HR–,
#         baseline HR–
#
# Notes:
#   • This script is designed for flat output (no nested subfolders).
#   • Scenarios are encoded directly into RDS filenames.
#
##########################################################################

suppressPackageStartupMessages({ library(fs) })

# ---------------------------------------------------------------------
# Construct root folder for this batch run: BC_Survivors_Batch_YYYYMMDD
# ---------------------------------------------------------------------
batch_root_dir <- function(date_obj = Sys.Date()) {
  file.path(getwd(), paste0("BC_Survivors_Batch_", format(date_obj, "%Y%m%d")))
}

# Numeric formatting helper for filenames
.frmt_num <- function(x, digits = 2) format(round(x, digits), nsmall = digits, trim = TRUE)

# ---------------------------------------------------------------------
# Filename generator encoding all scenario parameters
# ---------------------------------------------------------------------
scenario_filename <- function(bc_dx_start, bc_dx_end,
                              start_age_min, start_age_max,
                              end_age, early_sens,
                              OMST, DMST) {
  paste0(
    "diag", bc_dx_start, "-", bc_dx_end, "_",
    "start", start_age_min, "-", start_age_max, "_",
    "end", end_age, "_",
    "early", .frmt_num(early_sens), "_",
    "OMST", OMST, "_DMST", .frmt_num(DMST), ".rds"
  )
}

# =====================================================================
# run_and_save_flat()
#
# Description:
#   Runs ONE scenario (1 early_sens × 1 OMST/DMST pair)
#   across FOUR cohorts:
#     (a) HR+ Survivors
#     (b) HR+ Average-risk (baseline)
#     (c) HR– Survivors
#     (d) HR– Average-risk (baseline)
#
#   Saves all four results in a SINGLE .rds file.
#
# =====================================================================
run_and_save_flat <- function(
    bc_dx_start, bc_dx_end, years_to_mced,
    end_age,
    screen_int      = 1,
    early_sens,
    late_sens,
    OMST, DMST,
    fitted_data,
    qx_survivor_by_hr,
    qx_baseline,
    entry_weights_by_hr,
    cancer_list,
    date_obj = Sys.Date()
) {
  
  # -------------------------------------------------------------
  # 1. Compute MCED screening entry ages using the BC diagnosis window
  # -------------------------------------------------------------
  start_age_min <- bc_dx_start + years_to_mced
  start_age_max <- bc_dx_end   + years_to_mced
  start_age_seq <- start_age_min:start_age_max
  
  # Ensure output directory exists
  root_dir <- batch_root_dir(date_obj)
  dir_create(root_dir)
  
  # -------------------------------------------------------------
  # Helper to run MCED model for a specific HR cohort
  # -------------------------------------------------------------
  run_cell <- function(hr_status, cohort_type, qx_table, weights_df) {
    mced_multi_cancer_screening(
      hr_status     = hr_status,
      cohort_type   = cohort_type,
      start_age_seq = start_age_seq,
      end_age       = end_age,
      screen_int    = screen_int,
      early_sens    = early_sens,
      late_sens     = late_sens,
      fitted_data   = fitted_data,
      OMST          = OMST,
      DMST          = DMST,
      qx_table      = qx_table,
      cancer_list   = cancer_list,
      entry_weights = weights_df
    )
  }
  
  # Weight validation helper
  assert_weights <- function(df, tag) {
    if (!is.data.frame(df)) stop(sprintf("%s weights must be a data.frame", tag))
    if (!("age" %in% names(df))) stop(sprintf("%s weights must contain 'age'", tag))
    if (!any(c("w","weight") %in% names(df))) stop(sprintf("%s weights must contain 'w' or 'weight'", tag))
  }
  
  assert_weights(entry_weights_by_hr[["HR+"]], "HR+")
  assert_weights(entry_weights_by_hr[["HR-"]], "HR-")
  
  # -------------------------------------------------------------
  # 2. Run the four cohort × HR combinations
  # -------------------------------------------------------------
  res_surv_pos <- run_cell("HR+", "survivors", qx_survivor_by_hr[["HR+"]], entry_weights_by_hr[["HR+"]])
  res_base_pos <- run_cell("HR+", "baseline",  qx_baseline,                 entry_weights_by_hr[["HR+"]])
  res_surv_neg <- run_cell("HR-", "survivors", qx_survivor_by_hr[["HR-"]], entry_weights_by_hr[["HR-"]])
  res_base_neg <- run_cell("HR-", "baseline",  qx_baseline,                 entry_weights_by_hr[["HR-"]])
  
  # -------------------------------------------------------------
  # 3. Assemble payload
  # -------------------------------------------------------------
  payload <- list(
    meta = list(
      bc_dx_start    = bc_dx_start,
      bc_dx_end      = bc_dx_end,
      years_to_mced  = years_to_mced,
      start_age_min  = start_age_min,
      start_age_max  = start_age_max,
      end_age        = end_age,
      screen_int     = screen_int,
      early_sens     = early_sens,
      late_sens      = late_sens,
      OMST           = OMST,
      DMST           = DMST,
      cancer_list    = cancer_list
    ),
    survivors_HRpositive = res_surv_pos,
    baseline_HRpositive  = res_base_pos,
    survivors_HRnegative = res_surv_neg,
    baseline_HRnegative  = res_base_neg
  )
  
  # -------------------------------------------------------------
  # 4. Save the scenario as a flat RDS file
  # -------------------------------------------------------------
  out_file <- file.path(
    root_dir,
    scenario_filename(
      bc_dx_start, bc_dx_end,
      start_age_min, start_age_max,
      end_age, early_sens,
      OMST, DMST
    )
  )
  saveRDS(payload, out_file)
  
  invisible(out_file)
}

# =====================================================================
# run_all_scenarios()
#
# Description:
#   Outer batch runner:
#       • Loops over early_sens_vec
#       • Loops over OMST–DMST pairs
#       • Calls run_and_save_flat() for each combination
#       • Produces 6 flat files (if 3 pairs × 2 sens)
#
# =====================================================================
run_all_scenarios <- function(
    bc_dx_start, bc_dx_end, years_to_mced,
    end_age,
    screen_int      = 1,
    early_sens_vec,
    late_sens,
    OMST_DMST_pairs,
    fitted_data,
    qx_survivor_by_hr,
    qx_baseline,
    mced_start_age_dist_by_hr,
    cancer_list,
    date_obj = Sys.Date()
) {
  
  dir_create(batch_root_dir(date_obj))
  
  outputs <- c()
  
  for (es in early_sens_vec) {
    for (p in OMST_DMST_pairs) {
      
      outputs <- c(outputs,
                   run_and_save_flat(
                     bc_dx_start         = bc_dx_start,
                     bc_dx_end           = bc_dx_end,
                     years_to_mced       = years_to_mced,
                     end_age             = end_age,
                     screen_int          = screen_int,
                     early_sens          = es,
                     late_sens           = late_sens,
                     OMST                = p$OMST,
                     DMST                = p$DMST,
                     fitted_data         = fitted_data,
                     qx_survivor_by_hr   = qx_survivor_by_hr,
                     qx_baseline         = qx_baseline,
                     entry_weights_by_hr = mced_start_age_dist_by_hr,
                     cancer_list         = cancer_list,
                     date_obj            = date_obj
                   )
      )
      
    }
  }
  
  invisible(outputs)
}