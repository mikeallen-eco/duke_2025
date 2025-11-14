# analyze grassland bird density at Duke Farms 2010-2025 with Bayesian hierarchical distance sampling

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
                         mod = "bugs/distance_model_2025_pretty.txt")
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


### --- Step 4. get 2025 densities + change from previous year

(grsp_density_stats <- get_HDS_density(mod = "output/grsp_2025.rds"))
(bobo_density_stats <- get_HDS_density(mod = "output/bobo_2025.rds"))
(eame_density_stats <- get_HDS_density(mod = "output/eame_2025.rds"))
(rwbl_density_stats <- get_HDS_density(mod = "output/rwbl_2025.rds"))


# Make tables

grsp_table <- jags_table_html(readRDS("output/grsp_2025.rds"),full_table = F)
bobo_table <- jags_table_html(readRDS("output/bobo_2025.rds"))
eame_table <- jags_table_html(readRDS("output/eame_2025.rds"))
rwbl_table <- jags_table_html(readRDS("output/rwbl_2025.rds"))

kableExtra::save_kable(grsp_table, file = "grsp_jags_table.html")
kableExtra::save_kable(bobo_table, file = "bobo_jags_table.html")
kableExtra::save_kable(eame_table, file = "eame_jags_table.html")
kableExtra::save_kable(rwbl_table, file = "rwbl_jags_table.html")

### --- Step 6. plot density over time based on raw counts

# format raw count data

plot_raw_density 
ifelse()

raw_plot_data <- d %>% 
  filter(species %in% alpha)
  filter(!is.na(dist)) %>% 
  group_by(field, year, pt) %>% 
  tally() %>%
  right_join(allpts %>% select(year, field, pt) %>% distinct(), 
             by = join_by(pt, year, field)) %>%
  replace_na(list(n = 0)) %>%
  arrange(pt) %>%
  mutate(dens = n/(4*3.1415)) %>% 
  group_by(field, year) %>% 
  summarize(mean.num = mean(dens),
            sd = sd(dens),
            .groups = "drop") %>%
  mutate(npts = c(rep(10,6), rep(9,6)),
         yr = c(2010, 2012, 2013, 2018, 2019, 2024, 
                2010, 2012, 2013, 2018, 2019, 2024),
         se = sd/sqrt(npts))

# plot


raw.b <- ggplot(raw.plot) +
    geom_point(aes(x = yr, y = mean.num, color = field), size = 4,
               position = position_dodge(0.5)) +
    geom_errorbar(aes(x = yr, ymin = mean.num-(2*se), ymax = mean.num+(2*se),
                      color = field),
                  width = 0,
                  position = position_dodge(0.5)) +
    geom_line(aes(x = yr, y = mean.num, color = field),
              position = position_dodge(0.5)) +
    scale_color_manual(values = c("steelblue", "darkred")) +
    scale_x_continuous(breaks = seq(2010,2024, by = 4)) +
    ylim(0,2.5) +
    labs(y = "BOBO density (raw no./ha)", x = "", 
         color = "Field") +
    theme_bw() +
    theme(text = element_text(size = 13))
)

ggsave("figures/BOBO_6years.raw.field.year.png", height = 4, width = 6, dpi = 400)

