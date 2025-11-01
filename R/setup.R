# load libraries and functions

library(dplyr, quietly = T, warn.conflicts = F)
source("R/step1A_format_2010_2013_data_for_analysis.R")
source("R/step1B_format_2018_2019_data_for_analysis.R")
source("R/step1C_format_2024_2025_data_for_analysis.R")
source("R/step1D_get_jags_data_for_distance_sampling.R")
source("R/step2_run_JAGS_mod.R")
