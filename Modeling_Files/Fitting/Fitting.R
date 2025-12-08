#################################################################
#  Fitting.R
#  
#  Author: Kemal Caglar Gogebakan
#
#  This script fits age and stage-specific incidence to natural-history models 
#  for every cancer type
#  included in the multi-cancer early detection (MCED) screening
#  analysis (i.e., 16 cancer sites).
#
#  For each cancer site, and for each OMST–DMST pair, we:
#      • read SEER incidence data,
#      • apply a SIR multiplier (HR+ or HR− survivors),
#      • fit a natural-history model (single cancer),
#      • save metadata and fitted model objects in .Rdata files.
#
#  These fitted objects are later loaded by `load_fitted_data()`
#  and used by the MCED screening module to compute 
#  late-stage incidence reductions.
#
#  The functions below perform the following:
#
#    read_data():
#        Reads SEER incidence data
#
#    process_cancer_list():
#        Applies read_data() to each cancer site, using the
#        appropriate SIR multiplier for HR-positive,
#        HR-negative, or average-risk (SIR = 1).
#
#    fit_and_save_models():
#        For each cancer site:
#           - constructs metadata for all OMST–DMST pairs,
#           - calls get_fit() for each pair,
#           - stores the fitted results (<site>out_<label>)
#           - saves both metadata and results to disk.
#
#    get_fit():
#        Wrapper that:
#           - runs maximum-likelihood fitting (via get_max_LL),
#           - generates observed vs. expected plots,
#           - returns estimated rates, parameters, and summaries.
#
# =============================================================
library(here)
library(openxlsx)

# -------------------------------------------------------------
# Load functions used for fitting
# -------------------------------------------------------------
codepath <- here("Modeling_Files","Fitting")
source(file.path(codepath, "natural_history_code.R"))

# -------------------------------------------------------------
# Load cancer incidence data
# -------------------------------------------------------------
filepath <- here("Data", "MCED_Cancer_List_US_2011_2015.xlsx")
sheets <- openxlsx::getSheetNames(filepath)
cancer_list <- lapply(sheets, openxlsx::read.xlsx, xlsxFile = filepath)
names(cancer_list) <- sheets

# -------------------------------------------------------------
# Output folder: Modeling Files/Fitting
# -------------------------------------------------------------
filepath2 <- here("Modeling_Files", "Fitting", "Fitted_SEER_Data")

# Create dir if needed
if (!dir.exists(filepath2)) dir.create(filepath2, recursive = TRUE)

# -------------------------------------------------------------
# SIR multipliers for HR+ and HR− breast cancer survivors
# -------------------------------------------------------------
HR_pos_SIRs <- c(
  anus = 1.06, bladder = 1.07, breast = 1.71, cervix = 0.575, colorectal = 0.93,
  esophagus = 1.08, head = 1.37, kidney = 1.045, liver = 0.9, lung = 1.085,
  lymphoma = 0.965, melanoma = 1.15, ovary = 0.915, pancreas = 1.185,
  stomach = 1.15, thyroid = 1.2, uterus = 1.14
)

HR_neg_SIRs <- c(
  anus = 1.05, bladder = 0.935, breast = 2.07, cervix = 0.855, colorectal = 1.1,
  esophagus = 1.16, head = 1.63, kidney = 0.755, liver = 0.84, lung = 1.175,
  lymphoma = 0.965, melanoma = 0.995, ovary = 1.89, pancreas = 1.08,
  stomach = 1.19, thyroid = 1.08, uterus = 1.18
)

# -------------------------------------------------------------
# Read and process SEER data for one cancer site
# -------------------------------------------------------------
read_data <- function(the_data, multiplier = 1) {
  the_data$PY   <- the_data$Pop * the_data$years
  the_data$rate <- multiplier * the_data$Count / the_data$PY * 1e5
  the_data$Count <- the_data$Count * multiplier
  the_data <- the_data[the_data$midage > 0 & the_data$midage < 80, ]
  return(the_data)
}

# -------------------------------------------------------------
# Process entire cancer list with SIR multipliers
# -------------------------------------------------------------
process_cancer_list <- function(SIRs) {
  out_list <- list()
  for (site in names(cancer_list)) {
    site_clean <- tolower(gsub(" ", "_", site))
    dat <- cancer_list[[site]]
    multiplier <- ifelse(!is.na(SIRs[site_clean]), SIRs[site_clean], 1)
    out_list[[site_clean]] <- read_data(dat, multiplier)
  }
  return(out_list)
}

# -------------------------------------------------------------
# OMST–DMST combinations to fit
# -------------------------------------------------------------
OMST_DMST_pairs <- list(c(1, 0.5), c(2, 0.5), c(2, 1))
sites <- names(cancer_list)

# -------------------------------------------------------------
# Fit models for one group (HR+, HR−, or baseline)
# -------------------------------------------------------------
fit_and_save_models <- function(cancer_list_group, label) {
  
  for (site in sites) {
    site_clean <- tolower(gsub(" ", "_", site))
    
    OMST_site <- sapply(OMST_DMST_pairs, `[[`, 1)
    DMST_site <- sapply(OMST_DMST_pairs, `[[`, 2)
    
    metadata_var <- paste0("metadata_", label, "_", site_clean)
    assign(metadata_var, data.frame(
      fitID = seq_along(OMST_site),
      site  = site,
      OMST  = OMST_site,
      DMST  = DMST_site
    ), envir = .GlobalEnv)
    
    mean_value <- 3
    k_value    <- 16
    
    output_var <- paste0(site_clean, "out_", label)
    assign(output_var, mapply(
      FUN = "get_fit",
      OMST_site,
      DMST_site,
      MoreArgs = list(
        the_data   = cancer_list_group[[site_clean]],
        num_seeds  = 1,
        k          = k_value,
        mean1      = mean_value
      ),
      SIMPLIFY = FALSE
    ), envir = .GlobalEnv)
    
    save(
      list = c(output_var, metadata_var),
      file = file.path(filepath2, paste0(site_clean, "outfile_", label, ".Rdata"))
    )
  }
}

# -------------------------------------------------------------
# Fit all three groups (HR+, HR−, baseline)
# -------------------------------------------------------------
cancer_list_hr_pos  <- process_cancer_list(HR_pos_SIRs)
cancer_list_hr_neg  <- process_cancer_list(HR_neg_SIRs)
cancer_list_hr_base <- process_cancer_list(rep(1, length(cancer_list)))

fit_and_save_models(cancer_list_hr_pos,  label = "hr_pos")
fit_and_save_models(cancer_list_hr_neg,  label = "hr_neg")
fit_and_save_models(cancer_list_hr_base, label = "hr_base")