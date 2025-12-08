##############################################################################
# 01.data_read.R
#
# Author: Kemal Caglar Gogebakan
#
# Purpose:
#   Functions for loading and risk-adjusting breast cancer
#   survivor incidence and mortality data, as well as baseline population risk
#   inputs for MCED modeling.
#
# Description:
#   • Loads breast cancer (BC) non-metastatic incidence data (HR+ and HR−),
#     used to generate the diagnosis-age distribution for the BC survivor cohort.
#   • Loads female all-cause mortality rates.
#   • Loads all-cancers–combined incidence rates, used only for the first
#     five years after BC diagnosis (before MCED screening eligibility).
#   • Applies hazard ratios to adjust all-cause mortality for BC survivors.
#   • Applies SIR multipliers to inflate all-cancers–combined incidence for BC survivors.
#   • Loads baseline all-cause mortality for the average-risk population.
#   • Specifies the cancer types included in the MCED screening model.
#   • Defines data folders and file paths used by downstream functions.
#
# Contents:
#   .normalize_rate_table()  – Cleans and standardizes age/rate data.
#   .coerce_hr_vec()         – Ensures HR values are in a 2-element vector.
#   get_config()             – Central configuration (paths, SIRs, HRs, site list).
#   load_and_adjust_data()   – Main loader producing all model-ready inputs.
###################################################################################


suppressPackageStartupMessages({
  library(readxl)
  library(here)
})


# -------------------------------------------------------------------------
# Helper function: .normalize_rate_table
# -------------------------------------------------------------------------
# Converts any table containing an age column and a single-year rate column
# into a clean 2-column data frame (age, singleyears_rate), sorted by age.
.normalize_rate_table <- function(df, age_col = "age", rate_col = "singleyears_rate") {
  stopifnot(all(c(age_col, rate_col) %in% names(df)))
  out <- data.frame(
    age = as.integer(df[[age_col]]),
    singleyears_rate = as.numeric(df[[rate_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$age), ]
  rownames(out) <- NULL
  out
}


# -------------------------------------------------------------------------
# Helper function: .coerce_hr_vec
# -------------------------------------------------------------------------
# Converts HR input (scalar, list, or named vector) into a standardized,
# two-element numeric vector named HR+ and HR−.
.coerce_hr_vec <- function(x) {
  if (is.list(x)) x <- unlist(x)
  if (length(x) == 1L)
    return(c("HR+" = as.numeric(x), "HR-" = as.numeric(x)))
  stopifnot(all(c("HR+","HR-") %in% names(x)))
  as.numeric(x[c("HR+","HR-")]) |> setNames(c("HR+","HR-"))
}


# -------------------------------------------------------------------------
# Configuration: SIRs, HRs, file paths, cancer site list
# -------------------------------------------------------------------------
# Returns a list of file paths, hazard ratios, SIR values, data folder,
# and cancer sites to be modeled for MCED screening.
get_config <- function(
    SIR_0_5_pos = 1.12,
    SIR_0_5_neg = 1.38,
    HR_pos_acm  = 1.79,
    HR_neg_acm  = 1.79,
    data_dir    = here("Data"),
    fitted_dir  = here("Fitted_data"),
    cancer_list = c(
      "anus","bladder","cervix","colorectal","esophagus","head","kidney",
      "liver","lung","lymphoma","melanoma","ovary","pancreas","stomach",
      "thyroid","uterus"
    )
) {
  
  list(
    file_paths = list(
      bc_incidence_nonmet   = file.path(data_dir, "BC_incidence_nonmet_HR_2011_2015.xlsx"),
      all_cancers           = file.path(data_dir, "AllCancers_Incidence_2011_2015.xlsx"),
      all_cause_mort_female = file.path(data_dir, "Female_all_cause_mortality_2018.xlsx")
    ),
    SIRs = list(
      SIR_0_5_pos = SIR_0_5_pos,
      SIR_0_5_neg = SIR_0_5_neg
    ),
    HR_all_cause = c("HR+" = HR_pos_acm, "HR-" = HR_neg_acm),
    rdata_folder = fitted_dir,
    cancer_list  = cancer_list
  )
}


# -------------------------------------------------------------------------
# Main data loader and risk-adjuster function
# -------------------------------------------------------------------------
# Loads Excel files, standardizes them, applies risk adjustment for BC survivors,
# and returns:
#
#   • risk_adjusted:
#       - bc_nonmet_inc: incidence of non-metastatic BC (HR+ / HR−)
#       - qx_mortality_by_hr: all-cause mortality adjusted by HR
#       - all_cancers_inc_by_hr: any-cancer incidence adjusted by SIR
#
#   • baseline_risk:
#       - qx_mortality_by_hr: all-cause mortality for the average-risk population.
#         (Average-risk women do not have HR status; the HR+ and HR− entries are
#          identical and included only for easy comparison with survivor population.)

#   • rdata_folder: folder path for natural history model fitted data for 16 cancers for BC survivors and average-risk women
#   • cancer_list: vector of cancer sites to be included in MCED testing
#
load_and_adjust_data <- function(include_baseline = TRUE, cfg = get_config()) {
  fp   <- cfg$file_paths
  SIRs <- cfg$SIRs
  HRv  <- .coerce_hr_vec(cfg$HR_all_cause)
  
  # --- read input tables ----
  BC_nonmet_hr_pos_inc <- read_excel(fp$bc_incidence_nonmet,   sheet = "HR+")
  BC_nonmet_hr_neg_inc <- read_excel(fp$bc_incidence_nonmet,   sheet = "HR-")
  all_cancers_inc_raw  <- read_excel(fp$all_cancers,           sheet = "All_Cancers")
  female_acm_raw       <- read_excel(fp$all_cause_mort_female, sheet = "acm")
  
  # --- standardize to single-year tables ---
  all_cancers_inc <- .normalize_rate_table(all_cancers_inc_raw)
  female_acm      <- .normalize_rate_table(female_acm_raw)
  
  # --- survivor all-cause mortality: ACM × HR ---
  acm_pos <- transform(female_acm, singleyears_rate = singleyears_rate * HRv["HR+"])
  acm_neg <- transform(female_acm, singleyears_rate = singleyears_rate * HRv["HR-"])
  
  # --- survivor all-cancer incidence: incidence × SIR ---
  all_cancers_inc_hr_pos <- transform(
    all_cancers_inc, singleyears_rate = singleyears_rate * SIRs$SIR_0_5_pos
  )
  all_cancers_inc_hr_neg <- transform(
    all_cancers_inc, singleyears_rate = singleyears_rate * SIRs$SIR_0_5_neg
  )
  
  # --- packaged risk-adjusted bundle ---
  risk_adjusted <- list(
    bc_nonmet_inc         = list("HR+" = BC_nonmet_hr_pos_inc, "HR-" = BC_nonmet_hr_neg_inc),
    qx_mortality_by_hr    = list("HR+" = acm_pos,              "HR-" = acm_neg),
    all_cancers_inc_by_hr = list("HR+" = all_cancers_inc_hr_pos, "HR-" = all_cancers_inc_hr_neg)
  )
  
  # --- baseline population (optional) ---
  baseline_risk <- NULL
  if (include_baseline) {
    baseline_risk <- list(
      qx_mortality_by_hr = list("HR+" = female_acm, "HR-" = female_acm)
    )
  }
  
  # --- return everything ---
  list(
    risk_adjusted = risk_adjusted,
    baseline_risk = baseline_risk,
    rdata_folder  = cfg$rdata_folder,
    cancer_list   = cfg$cancer_list
  )
}