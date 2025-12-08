#######################################################################
# 06.multi_cancer_screening.R
#
# Author: Kemal Caglar Gogebakan
#
# Purpose:
#   Implement multi-cancer screening for:
#       • Breast cancer survivors (HR-positive / HR-negative)
#       • Age-matched average-risk women (“baseline”)
#
#   For each cancer site:
#       • Run the single-site MCED screening model
#       • Aggregate across start ages using entry-age weights
#       • Combine results across cancers into overall totals
#
# Functions:
# mced_multi_cancer_screening(..)
# Inputs:
#   hr_status       : "HR+" or "HR-"
#   cohort_type     : "survivors" or "baseline"
#                     (baseline uses HR=base fitted natural history models)
#   start_age_seq   : vector of start ages (e.g., 50:64)
#   end_age         : screening stops before this age
#   screen_int      : annual screening (must be 1)
#   early_sens      : MCED early-stage sensitivity (scalar)
#   late_sens       : MCED late-stage sensitivity (scalar)
#   fitted_data     : output from load_fitted_data()
#   OMST, DMST      : OMST–DMST pair selecting the fitted natural history model
#   qx_table        : all-cause mortality table (data.frame(age, singleyears_rate))
#   cancer_list     : ordered list of cancer sites
#   entry_weights   : data.frame(age, w) or (age, weight); weights sum to 1
#
# Output:
#   A list with:
#       $by_cancer_age   : list[cancer][start_age] -> 1-row df (unweighted)
#       $by_cancer       : weighted 1-row tibble per cancer (events, rates)
#       $overall         : overall totals summed across cancers
#
# Notes:
#   • baseline cohorts use “HR=base” natural history fits but retain
#     hr_status labels to preserve HR+ / HR– stratified comparisons.
##############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

# ---------------------------------------------------------------
# Multi-cancer screening function
# ---------------------------------------------------------------
mced_multi_cancer_screening <- function(
    hr_status,                          # "HR+" or "HR-"
    cohort_type = c("survivors","baseline"),
    start_age_seq,                      # e.g., 50:64
    end_age,                            # last screening age (exclusive)
    screen_int = 1,                     # annual screening only
    early_sens,                         # scalar
    late_sens,                          # scalar
    fitted_data,                        # from load_fitted_data()
    OMST, DMST,                         # chosen fitted pair
    qx_table,                           # ACM for HR group
    cancer_list,                        # cancer sites to evaluate
    entry_weights                       # df(age, w or weight)
) {
  
  cohort_type <- match.arg(cohort_type)
  
  # -------------------------------------------------------------
  # (0) Select model source:
  #     For survivors: HR+ or HR-
  #     For baseline: use HR=base natural history models
  # -------------------------------------------------------------
  model_hr_source <- if (cohort_type == "baseline") "HR=base" else hr_status
  fitted_for_run  <- fitted_data
  
  if (model_hr_source != hr_status) {
    # overwrite HR+ or HR– slot with HR=base fitted objects
    fitted_for_run$metadata_list[[hr_status]] <- fitted_data$metadata_list[[model_hr_source]]
    fitted_for_run$rateout_list [[hr_status]] <- fitted_data$rateout_list [[model_hr_source]]
  }
  
  # -------------------------------------------------------------
  # (1) Run single-cancer model for each cancer site
  # -------------------------------------------------------------
  by_cancer_age <- list()
  
  for (site in cancer_list) {
    
    # Optional progress printout
    if (isTRUE(getOption("mced.verbose", FALSE))) {
      ns_min <- max(0L, end_age - max(start_age_seq))
      ns_max <- max(0L, end_age - min(start_age_seq))
      message(sprintf("[MCED] %-14s | ages %d–%d | screens/start %d–%d",
                      site, min(start_age_seq), max(start_age_seq), ns_min, ns_max))
    }
    
    by_cancer_age[[site]] <- single_cancer_age_screening(
      cancer_site   = site,
      hr_status     = hr_status,
      start_age_seq = start_age_seq,
      end_age       = end_age,
      screen_int    = screen_int,
      early_sens    = early_sens,
      late_sens     = late_sens,
      fitted_data   = fitted_for_run,
      OMST1         = OMST,
      DMST1         = DMST,
      qx_table      = qx_table
    )
  }
  
  # -------------------------------------------------------------
  # (2) Process multi-cancer screening entry-age weights
  # -------------------------------------------------------------
  if (!is.data.frame(entry_weights))
    stop("entry_weights must be a data.frame")
  
  if (!("age" %in% names(entry_weights)))
    stop("entry_weights must contain column 'age'")
  
  w_col <- if ("w" %in% names(entry_weights)) {
    "w"
  } else if ("weight" %in% names(entry_weights)) {
    "weight"
  } else {
    stop("entry_weights must contain 'w' or 'weight'")
  }
  
  w_map <- setNames(as.numeric(entry_weights[[w_col]]),
                    as.integer(entry_weights$age))
  
  # align to start_age_seq
  w_map <- w_map[as.character(start_age_seq)]
  w_map[is.na(w_map)] <- 0
  
  # -------------------------------------------------------------
  # (3) Compute common weighted person-years denominator
  #     Using person_years for ANY cancer because PY depends
  #     only on start_age and ACM.
  # -------------------------------------------------------------
  first_site <- cancer_list[which(cancer_list %in% names(by_cancer_age))[1]]
  
  template_df <- bind_rows(by_cancer_age[[first_site]]) %>%
    mutate(w = w_map[as.character(start_age)])
  
  weighted_py <- sum(template_df$person_years * template_df$w)
  
  # -------------------------------------------------------------
  # (4) Weighted aggregation within each cancer across start ages
  # -------------------------------------------------------------
  by_cancer_tbl <- lapply(by_cancer_age, function(site_list) {
    
    df <- bind_rows(site_list) %>%
      mutate(w = w_map[as.character(start_age)])
    
    ctrl <- sum(df$control_late_events * df$w)
    scrn <- sum(df$screen_late_events  * df$w)
    avrt <- sum(df$averted_events      * df$w)
    
    tibble(
      control_late_events      = ctrl,
      screen_late_events       = scrn,
      averted_events           = avrt,
      person_years             = weighted_py,
      rate_ctrl_per_py         = ctrl / weighted_py,
      rate_screen_per_py       = scrn / weighted_py,
      rate_averted_per_py      = avrt / weighted_py,
      
      # event counts per 100k
      control_late_per100k     = ctrl * 1e5,
      screen_late_per100k      = scrn * 1e5,
      averted_per100k          = avrt * 1e5,
      
      # rates per 100k person-years
      rate_ctrl_per_100kPY     = (ctrl / weighted_py) * 1e5,
      rate_screen_per_100kPY   = (scrn / weighted_py) * 1e5,
      rate_averted_per_100kPY  = (avrt / weighted_py) * 1e5
    )
  })
  
  by_cancer <- by_cancer_tbl
  by_cancer_stacked <- bind_rows(by_cancer_tbl, .id = "cancer")
  
  # -------------------------------------------------------------
  # (5) Overall totals: sum numerators; keep common denominator
  # -------------------------------------------------------------
  total_ctrl <- sum(by_cancer_stacked$control_late_events)
  total_scrn <- sum(by_cancer_stacked$screen_late_events)
  total_avrt <- sum(by_cancer_stacked$averted_events)
  
  overall <- tibble(
    hr_status                 = hr_status,
    cohort_type               = cohort_type,
    control_late_events       = total_ctrl,
    screen_late_events        = total_scrn,
    averted_events            = total_avrt,
    person_years              = weighted_py,
    rate_ctrl_per_py          = total_ctrl / weighted_py,
    rate_screen_per_py        = total_scrn / weighted_py,
    rate_averted_per_py       = total_avrt / weighted_py,
    
    control_late_per100k      = total_ctrl * 1e5,
    screen_late_per100k       = total_scrn * 1e5,
    averted_per100k           = total_avrt * 1e5,
    
    rate_ctrl_per_100kPY      = (total_ctrl / weighted_py) * 1e5,
    rate_screen_per_100kPY    = (total_scrn / weighted_py) * 1e5,
    rate_averted_per_100kPY   = (total_avrt / weighted_py) * 1e5
  )
  
  # Return full object
  list(
    by_cancer_age = by_cancer_age,     # list[cancer][age] -> unweighted rows
    by_cancer     = by_cancer,         # 1-row tibble per cancer
    overall       = overall            # across all cancers combined
  )
}