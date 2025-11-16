# make output directory if needed
if(!dir.exists("output")){dir.create("output")}

# load libraries and functions

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)
source("R/step1A_format_2010_2013_data_for_analysis.R")
source("R/step1B_format_2018_2019_data_for_analysis.R")
source("R/step1C_format_2024_2025_data_for_analysis.R")
source("R/step1D_get_jags_data_for_distance_sampling.R")
source("R/step2_run_JAGS_mod.R")
source("R/step3_plot_HDS_density.R")
source("R/step4_get_HDS_density_function.R")
source("R/step5_jags_table_html_function.R")
source("R/step6_plot_raw_density_function.R")
source("R/step7_map_mean_counts_2025_function.R")

# vector of required CRAN packages
pkgs <- c(
  "lubridate", "tidyr", "jagsUI", "scales",
  "sf", "ggplot2", "ggspatial", "patchwork",
  "knitr", "kableExtra"
)

# install any that are missing
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Package '%s' not found. Installing...", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(pkgs, install_if_missing))
