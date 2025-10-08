# subset and format 2018-2019 data for comparison with 2024-2025 points
# selects the subset of points that correspond spatially to the 2024-2025 points
# and formats them for distance sampling

library(dplyr, quietly = T, warn.conflicts = F)
library(lubridate, quietly = T, warn.conflicts = F)

format_2024_2025_data_for_analysis <- function(path = "data/Grassland_Birds_2024_data.csv"){
  
  # read in the 2024 or 2025 data (same input format)
  d <- read.csv(path) %>%
    dplyr::rename(observer = Collector) %>%
    mutate(year = paste0("Y", substr(year,3,4)),
           date = lubridate::mdy(date),
           time = hour(hm(start)) + minute(hm(start))/60,
           pt = paste0(year, pt),
           obsperpt = paste0(observer, period, pt),
           perpt = paste0(period, pt)) %>%
    select(observer, year, date, time, pt, field, period, obsperpt, perpt,
           species, num, dist = distance)
  
  return(d)
  
}