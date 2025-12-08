##########################################################################
# 05. single_cancer_screening.R
#
# Author: Kemal Ca
#
# Purpose:
#   Projects late-stage cancer outcomes for a *single cancer site*
#       - HR-positive breast cancer survivors
#       - HR-negative breast cancer survivors
#       - Average-risk women
#   using calibrated data by OMST, LMST and screening program design.
#
#   For each start age (e.g., 50–64):
#       • screening occurs annually until `end_age`
#       • individuals continue to be followed until death
#         (via external all-cause mortality) or age 100 (lifetime)
#       • stage-shift model projects late-stage cancer incidence
#         under screening and no screening
#
# Key tasks in this file:
#   • Compute person-years.
#   • Use all-cause mortality as a competing risk to cancer diagnosis during screening period
#   • Project late-stage incidence and rates per person-year and
#     relative reductions.
#
# Functions:
#   1) .compute_py(qx_table, start_age, end_age)
#       - Computes person-years.
#
#   2) single_cancer_age_screening(...)
#       - Main interface:
#           - Loops over start ages
#           - Selects fitted natural history model
#           - Projects screening vs control late-stage outcomes
#           - Diagnosis competes with all-cause mortality
#           - Returns a list of data frames (one per start age)
#
# Output structure:
#   A named list with one data frame per start age:
#       $`50`  -> results for start_age=50
#       $`51`
#       ...
#
#   Each data frame contains:
#       cancer_site
#       hr_status
#       start_age
#       early_sens
#       late_sens
#       control_late_events
#       screen_late_events
#       averted_events
#       person_years
#       rate_ctrl_per_py
#       rate_screen_per_py
#       rate_averted_per_py
#       relative_reduction_late_stage
#
####################################################################

suppressPackageStartupMessages({ library(dplyr) })

# ---------------------------------------------------------------
# Compute discrete person-years using survival S(a)
# ---------------------------------------------------------------
# Inputs:
#   qx_table  : data.frame(age, singleyears_rate)
#   start_age : starting age
#   end_age   : ending age (non-inclusive)
#
# Output:
#   Numeric value: sum_{a=start_age}^{end_age-1} S(a)
#
.compute_py <- function(qx_table, start_age, end_age) {
  if (start_age >= end_age) return(0)
  qx <- setNames(qx_table$singleyears_rate, as.character(qx_table$age))
  S  <- 1.0
  py <- 0.0
  
  for (a in seq.int(start_age, end_age - 1L)) {
    py <- py + S
    r <- qx[as.character(a)]
    if (is.na(r)) r <- 0
    r <- min(max(r, 0), 0.9999)
    S <- S * (1 - r)
  }
  py
}

# ---------------------------------------------------------------
# Main function: outcome evaluation for a single cancer site
# ---------------------------------------------------------------
# Inputs:
#   cancer_site   : string (e.g., "lung", "pancreas")
#   hr_status     : "HR+" or "HR-"
#   start_age_seq : vector of start ages (e.g., 50:64)
#   end_age       : screening stops at this age
#   screen_int    : assumed annual (=1)
#   early_sens    : early-stage sensitivity
#   late_sens     : late-stage sensitivity
#   fitted_data   : output from load_fitted_data()
#   OMST1, DMST1  : choose the specific fitted model
#   qx_table      : ACM rates (data.frame(age, singleyears_rate))
#
# Output:
#   Named list of data frames, one per start age.
# ---------------------------------------------------------------
single_cancer_age_screening <- function(
    cancer_site,
    hr_status,               # "HR+" or "HR-"
    start_age_seq,           # e.g., 50:64
    end_age,                 # screening stops at this age
    screen_int = 1,          # annual
    early_sens,              # scalar
    late_sens,               # scalar
    fitted_data,             # from load_fitted_data()
    OMST1, DMST1,            # select the fitted pair explicitly
    qx_table                 # data.frame(age, singleyears_rate)
) {
  stopifnot(screen_int == 1)  # only annual supported
  
  follow_end_age <- 100L      # follow-up until death (via ACM)
  
  # ---- Retrieve fitted natural history model ----
  meta    <- fitted_data$metadata_list[[hr_status]][[cancer_site]]
  ro_bank <- fitted_data$rateout_list  [[hr_status]][[cancer_site]]
  
  if (is.null(meta) || is.null(ro_bank)) {
    stop(sprintf("Missing fitted objects for %s / %s", cancer_site, hr_status))
  }
  
  sel <- meta[meta$OMST == OMST1 & meta$DMST == DMST1, ]
  if (nrow(sel) != 1L) {
    stop(sprintf("Could not uniquely select fit for %s / %s with OMST=%s, DMST=%s",
                 cancer_site, hr_status, as.character(OMST1), as.character(DMST1)))
  }
  
  rateout <- ro_bank[[ sel$fitID ]]
  
  # ---- Loop over start ages ----
  out <- vector("list", length(start_age_seq))
  names(out) <- as.character(start_age_seq)
  
  # Pre-index ACM for speed
  qx_vec_all <- setNames(qx_table$singleyears_rate, as.character(qx_table$age))
  
  for (start_age in start_age_seq) {
    
    # total lifetime horizon
    total_intervals <- max(0L, follow_end_age - start_age)
    
    # screening horizon
    screen_stop_age <- min(end_age, follow_end_age)
    n_screens <- max(0L, screen_stop_age - start_age)
    
    # Degenerate case: start age >= 100
    if (total_intervals == 0L) {
      py <- .compute_py(qx_table, start_age, follow_end_age)
      df <- data.frame(
        cancer_site  = cancer_site,
        hr_status    = hr_status,
        start_age    = start_age,
        early_sens   = early_sens,
        late_sens    = late_sens,
        control_late_events = 0,
        screen_late_events  = 0,
        averted_events      = 0,
        person_years        = py,
        rate_ctrl_per_py    = if (py > 0) 0 else NA_real_,
        rate_screen_per_py  = if (py > 0) 0 else NA_real_,
        rate_averted_per_py = if (py > 0) 0 else NA_real_,
        relative_reduction_late_stage = NA_real_,
        stringsAsFactors = FALSE
      )
      out[[as.character(start_age)]] <- df
      next
    }
    
    # No screening at all for this start_age
    if (n_screens == 0L) {
      py <- .compute_py(qx_table, start_age, follow_end_age)
      df <- data.frame(
        cancer_site  = cancer_site,
        hr_status    = hr_status,
        start_age    = start_age,
        early_sens   = early_sens,
        late_sens    = late_sens,
        control_late_events = 0,
        screen_late_events  = 0,
        averted_events      = 0,
        person_years        = py,
        rate_ctrl_per_py    = if (py > 0) 0 else NA_real_,
        rate_screen_per_py  = if (py > 0) 0 else NA_real_,
        rate_averted_per_py = if (py > 0) 0 else NA_real_,
        relative_reduction_late_stage = NA_real_,
        stringsAsFactors = FALSE
      )
      out[[as.character(start_age)]] <- df
      next
    }
    
    # Number of post-screen follow-up intervals
    n_followup <- total_intervals - n_screens
    stopifnot(n_followup >= 0L)
    
    # Run full stage-shift model up to lifetime horizon
    sh <- stage_shift_by_screen_single(
      numscreens             = n_screens,
      age                    = start_age,
      screen.int             = 1,
      rateout                = rateout,
      sens_e                 = early_sens,
      sens_l                 = late_sens,
      num.followup.intervals = n_followup
    )
    
    # Extract increments from cumulative curves
    inc_control   <- c(sh$cum_control_late[1],  diff(sh$cum_control_late))
    inc_screen    <- c(sh$cum_screen_late[1],   diff(sh$cum_screen_late))
    inc_interval  <- c(sh$cum_interval_late[1], diff(sh$cum_interval_late))
    inc_screen_total <- inc_screen + inc_interval
    
    stopifnot(length(inc_control) == total_intervals - 1L)
    
    # External ACM survival weighting
    ages_seq <- as.character(seq.int(start_age, follow_end_age - 1L))
    qx_seq   <- qx_vec_all[ages_seq]
    qx_seq[is.na(qx_seq)] <- 0
    qx_seq   <- pmin(pmax(qx_seq, 0), 0.9999)
    
    # S(a): survival to interval start
    S_to <- cumprod(c(1, 1 - qx_seq))[seq_len(total_intervals)]
    
    # Expected events
    E_ctrl <- sum(S_to * inc_control)
    E_scr  <- sum(S_to * inc_screen_total)
    PY     <- sum(S_to)
    
    # Output row
    df <- data.frame(
      cancer_site                 = cancer_site,
      hr_status                   = hr_status,
      start_age                   = start_age,
      early_sens                  = early_sens,
      late_sens                   = late_sens,
      control_late_events         = E_ctrl,
      screen_late_events          = E_scr,
      averted_events              = E_ctrl - E_scr,
      person_years                = PY,
      rate_ctrl_per_py            = if (PY > 0) E_ctrl / PY else NA_real_,
      rate_screen_per_py          = if (PY > 0) E_scr  / PY else NA_real_,
      rate_averted_per_py         = if (PY > 0) (E_ctrl - E_scr) / PY else NA_real_,
      relative_reduction_late_stage =
        if (E_ctrl > 0) (E_ctrl - E_scr) / E_ctrl else NA_real_,
      stringsAsFactors = FALSE
    )
    
    out[[as.character(start_age)]] <- df
  }
  
  out
}