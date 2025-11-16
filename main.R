# analyze grassland bird density at Duke Farms 2010-2025 with Bayesian hierarchical distance sampling

# Initialize (run this to load functions and libraries)
source("R/setup.R")


### --- Step 1. Read in and format data

# read in data for each set of years
d10_13 <- format_2010_2013_data_for_analysis("data/Grassland_Birds_2010-2013_data.csv")
d18_19 <- format_2018_2019_data_for_analysis("data/Grassland_Birds_2018-2019_data.csv")
d24 <- format_2024_2025_data_for_analysis("data/Grassland_Birds_2024_data.csv")
d25 <- format_2024_2025_data_for_analysis("data/20250520_grassland_birds_2025_Allen/Grassland_Birds_2025_data.csv")
d <- bind_rows(d10_13, d18_19, d24, d25)

# format JAGS data list for each species
b <- get_jags_data_for_distance_sampling(spname = "BOBO", dist_threshold = 100)
g <- get_jags_data_for_distance_sampling(spname = "GRSP", dist_threshold = 100)
e <- get_jags_data_for_distance_sampling(spname = "EAME", dist_threshold = 100)
r <- get_jags_data_for_distance_sampling(spname = "RWBL", dist_threshold = 100)


### --- Step 2. run JAGS models

# note: view and edit JAGS code file in: JAGS/distance_model_2025.txt

# --- Bobolink
set.seed(206)
bobo_out <- run_JAGS_mod(jags_data = b, 
                         sp_name = "bobo",
                         mod = "bugs/distance_model_2025.txt")
bobo_out # examine output
# jagsUI:: traceplot(bobo_out) # examine trace plots

# --- Grasshopper Sparrow
set.seed(206)
grsp_out <- run_JAGS_mod(jags_data = g, sp_name = "grsp",
                         mod = "bugs/distance_model_2025.txt")
grsp_out # examine output
# jagsUI:: traceplot(grsp_out) # examine trace plots

# --- Eastern Meadowlark
set.seed(206)
eame_out <- run_JAGS_mod(jags_data = e, sp_name = "eame",
                         mod = "bugs/distance_model_2025.txt")
eame_out # examine output
# jagsUI:: traceplot(eame_out) # examine trace plots

# --- Red-winged Blackbird
set.seed(206)
rwbl_out <- run_JAGS_mod(jags_data = r, sp_name = "rwbl",
                         mod = "bugs/distance_model_2025.txt")
rwbl_out # examine output
# jagsUI:: traceplot(eame_out) # examine trace plots


### --- Step 3. plot density over time based on hierarchical distance sampling model

(bobo_density <- plot_HDS_density(mod = "output/bobo_2025.rds"))
(grsp_density <- plot_HDS_density(mod = "output/grsp_2025.rds"))
(eame_density <- plot_HDS_density(mod = "output/eame_2025.rds"))
(rwbl_density <- plot_HDS_density(mod = "output/rwbl_2025.rds"))


### --- Step 4. get 2025 summary stats (counts, densities, & change from previous year)

# counts
(counts <- d %>%
  filter(year %in% "Y25",
         dist < 101,
         species %in% c("EAME", "BOBO", "GRSP", "RWBL")) %>%
  group_by(species, period, field) %>%
  summarize(sum = sum(num, na.rm = T)) %>%
    arrange(species, field, period))

(grsp_density_stats <- get_HDS_density(mod = "output/grsp_2025.rds"))
(bobo_density_stats <- get_HDS_density(mod = "output/bobo_2025.rds"))
(eame_density_stats <- get_HDS_density(mod = "output/eame_2025.rds"))
(rwbl_density_stats <- get_HDS_density(mod = "output/rwbl_2025.rds"))


# --- Step 5 - make tables for inclusion in the report

grsp_table <- jags_table_html(readRDS("output/grsp_2025.rds"),full_table = F)
bobo_table <- jags_table_html(readRDS("output/bobo_2025.rds"))
eame_table <- jags_table_html(readRDS("output/eame_2025.rds"))
rwbl_table <- jags_table_html(readRDS("output/rwbl_2025.rds"))

kableExtra::save_kable(grsp_table, file = "output/grsp_jags_table.html")
kableExtra::save_kable(bobo_table, file = "output/bobo_jags_table.html")
kableExtra::save_kable(eame_table, file = "output/eame_jags_table.html")
kableExtra::save_kable(rwbl_table, file = "output/rwbl_jags_table.html")


### --- Step 6. plot density over time based on raw counts

plot_raw_density(data = d, alpha = "GRSP")
plot_raw_density(data = d, alpha = "BOBO")
plot_raw_density(data = d, alpha = "EAME")
plot_raw_density(data = d, alpha = "RWBL")


### --- Step 7. Mapping mean counts

(grsp_map <- map_mean_counts(data = d, alpha = "GRSP"))
(bobo_map <- map_mean_counts(data = d, alpha = "BOBO"))
(eame_map <- map_mean_counts(data = d, alpha = "EAME"))
(rwbl_map <- map_mean_counts(data = d, alpha = "RWBL"))
library(patchwork)
(grsp_map | bobo_map) / (eame_map | rwbl_map)

 ggsave("output/all_map.png", height = 8, width = 12, dpi = 600)
