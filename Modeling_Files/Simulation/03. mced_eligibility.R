###################################################################
# 03. mced eligibility.R
#
# Author: Kemal Caglar Gogebakan
#
# Purpose:
#   Construct MCED-screening entry-age distributions for:
#     • Breast cancer survivors (HR-positive and HR-negative),
#       followed k years from diagnosis.
#     • Age-matched average-risk women, whose entry-age distributions
#       are directly copied from survivors.
#
# Conceptual flow:
#   1) Use non-metastatic breast cancer incidence data (HR-positive / HR-negative)
#      to create an age distribution at diagnosis (weights by age).
#   2) Follow each diagnosis-age cohort forward for Y years:
#        – Remove women who die of any cause using age-specific all-cause mortality
#          (already hazard-ratio–adjusted for survivors).
#        – Remove subsequent primary cancers using all-cancers incidence
#          (SIR-adjusted for survivors).
#   3) Convert the remaining eligible proportions into MCED screening entry-age
#      weight distributions (weights sum to 1 within each HR group).
#   4) Copy these survivor-based entry-age distributions to define
#      the age-matched average-risk MCED entry cohorts.
#
# Functions in this file:
#   • .to_prob(x)
#       Convert a rate per 100,000 into a probability if needed.
#
#   • .vals_by_age(obj, ages)
#       Extract age-indexed values from:
#         - data.frame(age, singleyears_rate), or
#         - named vector whose names are ages.
#
#   • .get_dx_age_distribution(df, diag_age_min, diag_age_max)
#       Derive normalized age distribution at breast cancer diagnosis
#       using counts, Rate×Pop, or singleyears_rate as the underlying input.
#
#   • .follow_k_years(ages0, init_weights, lambda_all_cancers, acm, k)
#       Follow each diagnosis-age cohort forward k years by applying:
#         Survival = Π (1 − ACM) × (1 − all-cancer incidence),
#       and return the normalized surviving age-distribution weights.
#
#   • build_mced_entry_ages(risk_adjusted, diag_age_min, diag_age_max, years_to_screen_start)
#       Main function that returns MCED entry-age distributions for
#       breast cancer survivors and age-matched baseline women by HR group.
#
# Output:
#   list(
#     survivors_entry_age = list(
#       "HR+" = data.frame(age, weight),
#       "HR-" = data.frame(age, weight)
#     ),
#     baseline_entry_age  = list(
#       "HR+" = data.frame(age, weight),
#       "HR-" = data.frame(age, weight)
#     )
#   )
#
# Notes:
#   - Average-risk women do not have hormone-receptor status; we retain
#     the HR+/HR− labels only to mirror the survivor structure and enable
#     age-matched comparisons in downstream modeling.
#   - Within each HR group, weight represents the normalized probability
#     of entering MCED screening at a given age.
# ##########################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

# ------------------------------------------------------------
# Convert rate per 100k into a probability if needed
# ------------------------------------------------------------
.to_prob <- function(x) {
  if (is.null(x)) return(NULL)
  if (any(x > 1, na.rm = TRUE)) x / 1e5 else x
}

# ------------------------------------------------------------
# Extract age-indexed values from df(age, singleyears_rate)
# or from a named vector
# ------------------------------------------------------------
.vals_by_age <- function(obj, ages) {
  if (is.null(obj)) return(rep(0, length(ages)))
  
  if (is.vector(obj) && !is.null(names(obj))) {
    return(as.numeric(obj[as.character(ages)]))
  }
  
  if (is.data.frame(obj)) {
    return(as.numeric(obj$singleyears_rate[match(ages, obj$age)]))
  }
  
  stop("Input must be a named vector by age or a data.frame(age, singleyears_rate).")
}

# ------------------------------------------------------------
# Derive diagnosis-age distribution (weights by age)
# from BC non-met incidence input.
#
# Accepts:
#   - Count
#   - Rate + Pop  → Rate × Pop
#   - singleyears_rate → used as proxy count if needed
# ------------------------------------------------------------
.get_dx_age_distribution <- function(df, diag_age_min, diag_age_max) {
  
  sub <- df[df$age >= diag_age_min & df$age <= diag_age_max, , drop = FALSE]
  sub <- sub[order(sub$age), ]
  if (!nrow(sub)) stop("No rows in the diagnosis-age window.")
  
  if ("Count" %in% names(sub)) {
    counts <- sub$Count
    
  } else if (all(c("Rate", "Pop") %in% names(sub))) {
    r <- if (any(sub$Rate > 1, na.rm = TRUE)) sub$Rate / 1e5 else sub$Rate
    counts <- r * sub$Pop
    
  } else if ("singleyears_rate" %in% names(sub)) {
    r <- if (any(sub$singleyears_rate > 1, na.rm = TRUE)) {
      sub$singleyears_rate / 1e5
    } else sub$singleyears_rate
    counts <- r
    
  } else {
    stop("Provide Count, or Rate+Pop, or singleyears_rate in the BC incidence table.")
  }
  
  # Normalize counts → weights
  w <- if (sum(counts) > 0) counts / sum(counts) else rep(0, length(counts))
  
  tibble(age = as.integer(sub$age), weight = as.numeric(w))
}

# ------------------------------------------------------------
# Follow k years forward:
#   Survival_i = Π (1 − ACM(age+y)) × (1 − IncAll(age+y))
#
# Inputs:
#   ages0         — diagnosis ages
#   init_weights  — normalized weights at diagnosis
#   lambda_all_cancers — all-cancer incidence df or vector
#   acm           — all-cause mortality df or vector
#   k             — years until MCED screening eligibility
#
# Returns: vector of normalized surviving weights
# ------------------------------------------------------------
.follow_k_years <- function(ages0, init_weights, lambda_all_cancers, acm, k) {
  
  surv <- sapply(seq_along(ages0), function(i) {
    
    yrs <- ages0[i]:(ages0[i] + k - 1)
    
    qx <- pmin(pmax(.vals_by_age(acm, yrs), 0), 1)  # ACM
    lam <- pmin(pmax(.to_prob(.vals_by_age(lambda_all_cancers, yrs)), 0), 1)  # All-cancers
    
    prod((1 - qx) * (1 - lam))
  })
  
  out <- init_weights * surv
  if (sum(out) > 0) out / sum(out) else out
}

# ------------------------------------------------------------
# Main interface
# ------------------------------------------------------------
build_mced_entry_ages <- function(risk_adjusted,
                                  diag_age_min,
                                  diag_age_max,
                                  years_to_screen_start = 5) {
  
  stopifnot(
    is.list(risk_adjusted),
    all(c("bc_nonmet_inc", "qx_mortality_by_hr", "all_cancers_inc_by_hr") %in% names(risk_adjusted))
  )
  
  # ------------------------
  # HR-positive survivors
  # ------------------------
  dx_pos <- .get_dx_age_distribution(
    df = risk_adjusted$bc_nonmet_inc[["HR+"]],
    diag_age_min = diag_age_min,
    diag_age_max = diag_age_max
  )
  
  w_pos_post <- .follow_k_years(
    ages0 = dx_pos$age,
    init_weights = dx_pos$weight,
    lambda_all_cancers = risk_adjusted$all_cancers_inc_by_hr[["HR+"]],
    acm = risk_adjusted$qx_mortality_by_hr[["HR+"]],
    k = years_to_screen_start
  )
  
  hr_pos_df <- tibble(
    age = dx_pos$age + years_to_screen_start,
    weight = w_pos_post
  )
  
  # ------------------------
  # HR-negative survivors
  # ------------------------
  dx_neg <- .get_dx_age_distribution(
    df = risk_adjusted$bc_nonmet_inc[["HR-"]],
    diag_age_min = diag_age_min,
    diag_age_max = diag_age_max
  )
  
  w_neg_post <- .follow_k_years(
    ages0 = dx_neg$age,
    init_weights = dx_neg$weight,
    lambda_all_cancers = risk_adjusted$all_cancers_inc_by_hr[["HR-"]],
    acm = risk_adjusted$qx_mortality_by_hr[["HR-"]],
    k = years_to_screen_start
  )
  
  hr_neg_df <- tibble(
    age = dx_neg$age + years_to_screen_start,
    weight = w_neg_post
  )
  
  # ------------------------
  # Baseline (average-risk)
  # Just copy survivors’ entry distributions (age-matched)
  # ------------------------
  list(
    survivors_entry_age = list("HR+" = hr_pos_df, "HR-" = hr_neg_df),
    baseline_entry_age  = list("HR+" = hr_pos_df, "HR-" = hr_neg_df)
  )
}