###########################################################################
# 02. load fitted data.R
#
# Author: Kemal Caglar Gogebakan
#
# Loads calibrated natural-history fit objects (metadata + rate outputs)
# for MCED screening across 16 cancer sites.
#
# These fitted objects correspond to the calibrated incidence models
# used in the MCED simulations for:
#   • Breast cancer survivors stratified by hormone-receptor status (HR+ / HR−)
#   • Age-matched average-risk women (HR=base)
#
# Each fitted object reflects the calibration performed for the OMST/DMST
# (mean sojourn time) parameter pairs used in the model.
#
# Expected Rdata filenames per cancer site:
#     <site>outfile_hr_pos.Rdata      # HR+ survivors
#     <site>outfile_hr_neg.Rdata      # HR− survivors
#     <site>outfile_hr_base.Rdata     # Average-risk
#
# Each file must contain two objects:
#     metadata_hr_<pos|neg|base>_<site>
#     <site>out_hr_<pos|neg|base>
#
# Returns:
#     list(
#       metadata_list = list(HR+ = ..., HR− = ..., HR=base = ...),
#       rateout_list  = list(HR+ = ..., HR− = ..., HR=base = ...)
#     )
###########################################################################

suppressPackageStartupMessages({ library(stringr) })

library(stringr)

load_fitted_data <- function(rdata_folder, cancer_list, verbose = TRUE) {
  hr_groups <- c("HR+" = "hr_pos", "HR-" = "hr_neg", "HR=base" = "hr_base")
  
  rdata_files <- list.files(path = rdata_folder, pattern = "\\.Rdata$", full.names = TRUE)
  
  metadata_list <- list("HR+" = list(), "HR-" = list(), "HR=base" = list())
  rateout_list  <- list("HR+" = list(), "HR-" = list(), "HR=base" = list())
  
  for (file in rdata_files) {
    file_base <- basename(file)
    
    # Match pattern like pancreasoutfile_hr_neg.Rdata or pancreasoutfile_hr_base.Rdata
    matches <- str_match(file_base, "^([a-z]+)outfile_hr_(neg|pos|base)\\.Rdata$")
    if (any(is.na(matches))) {
      warning(paste("Skipping unrecognized file:", file_base))
      next
    }
    
    cancer <- matches[2]
    if (!(cancer %in% cancer_list)) {
      next
    }
    
    hr_code <- paste0("hr_", matches[3])
    hr_status <- names(hr_groups)[hr_groups == hr_code]
    
    temp_env <- new.env()
    load(file, envir = temp_env)
    obj_names <- ls(temp_env)
    
    metadata_name <- paste0("metadata_", hr_code, "_", cancer)
    rateout_name  <- paste0(cancer, "out_", hr_code)
    
    if (!(metadata_name %in% obj_names)) {
      warning(paste("Missing metadata object in file:", file_base))
      next
    }
    if (!(rateout_name %in% obj_names)) {
      warning(paste("Missing rateout object in file:", file_base))
      next
    }
    
    metadata_list[[hr_status]][[cancer]] <- get(metadata_name, envir = temp_env)
    rateout_list[[hr_status]][[cancer]]  <- get(rateout_name,  envir = temp_env)
  }
  
  if (verbose) {
    for (hr in names(hr_groups)) {
      for (cancer in cancer_list) {
        if (is.null(metadata_list[[hr]][[cancer]])) {
          warning(paste("Missing metadata for", hr, cancer))
        }
        if (is.null(rateout_list[[hr]][[cancer]])) {
          warning(paste("Missing rateout for", hr, cancer))
        }
      }
    }
  }
  
  return(list(metadata_list = metadata_list, rateout_list = rateout_list))
}
