# analyze grassland bird density at Duke Farms 2010-2025 with hierarchical distance sampling

source("R/setup.R")

# --- Read in data

d10 <- format_2010_2013_data_for_analysis("data/Grassland_Birds_2010-2013_data.csv")
d18 <- format_2018_2019_data_for_analysis("data/Grassland_Birds_2018-2019_data.csv")
d24 <- format_2024_2025_data_for_analysis("data/Grassland_Birds_2024_data.csv")
d25 <- format_2024_2025_data_for_analysis("data/20250520_grassland_birds_2025_Allen/Grassland_Birds_2025_data.csv")

d <- bind_rows(d18,d24,d25)

survey_info <- d %>%
  select(period, pt, observer, field, year, date, time) %>%
  distinct() %>%
  arrange(pt, year, period) %>%
  mutate(perpt = paste0(period, pt)) %>%
  mutate(dup = duplicated(perpt)*1)
