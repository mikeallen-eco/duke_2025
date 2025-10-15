# analyze grassland bird density at Duke Farms 2010-2025 with Bayesian hierarchical distance sampling

source("R/setup.R")

### --- Step 1. Read in and format data

# read in data for each set of years
d10_13 <- format_2010_2013_data_for_analysis("data/Grassland_Birds_2010-2013_data.csv")
d18_19 <- format_2018_2019_data_for_analysis("data/Grassland_Birds_2018-2019_data.csv")
d24 <- format_2024_2025_data_for_analysis("data/Grassland_Birds_2024_data.csv")
d25 <- format_2024_2025_data_for_analysis("data/20250520_grassland_birds_2025_Allen/Grassland_Birds_2025_data.csv")
d <- bind_rows(d10_13, d18_19,d24,d25)

# format JAGS data list for each species
b <- get_jags_data_for_distance_sampling(spname = "BOBO", dist_threshold = 100)

### --- Step 2. Edit & save JAGS model

# JAGS model specification for point transect data
cat("
model {

  #### ------------------------
  #### Abundance model
  #### ------------------------

  # Priors
  beta_int       ~ dnorm(0, 0.04)      # sd ≈ 5
  beta_field     ~ dnorm(0, 0.25)      # sd = 2
  
  # Year effects
  beta_Y10       ~ dnorm(0, 0.25)
  beta_Y12       ~ dnorm(0, 0.25)
  beta_Y13       ~ dnorm(0, 0.25)
  beta_Y19       ~ dnorm(0, 0.25)
  beta_Y24       ~ dnorm(0, 0.25)
  beta_Y25       ~ dnorm(0, 0.25)
  
  # Field × Year interactions
  beta_fieldY10  ~ dnorm(0, 0.25)
  beta_fieldY12  ~ dnorm(0, 0.25)
  beta_fieldY13  ~ dnorm(0, 0.25)
  beta_fieldY19  ~ dnorm(0, 0.25)
  beta_fieldY24  ~ dnorm(0, 0.25)
  beta_fieldY25  ~ dnorm(0, 0.25)

  #### ------------------------
  #### Detection model
  #### ------------------------

  # Hyperpriors for yearly detection variance
  mu_alpha    ~ dnorm(3.5, 0.25)       # mean log(sigma) ≈ 3.5 (≈33 m)
  sigma_alpha ~ dunif(0, 3)
  tau_alpha   <- 1 / (sigma_alpha * sigma_alpha)

  # Year-level random effects
  for (y in 1:nyears) {
  alpha_year[y] ~ dnorm(mu_alpha, tau_alpha)
  alpha_year_used[y] <- ifelse(has_det[y] == 1, alpha_year[y], mu_alpha)
  }


  #### ------------------------
  #### Observation model (marginalized)
  #### ------------------------

  for (i in 1:nind) {
    dclass[i] ~ dcat(fc[site[i],])    # Part 1 of HM
  }

  for (s in 1:nsites) {

    # Year-specific sigma
    log(sigma[s]) <- mu_alpha + alpha_year_used[ year[s] ]

    # Construct distance band detection probabilities
    for (g in 1:nD) {
      log(p[s,g]) <- -midpt[g] * midpt[g] / (2 * sigma[s] * sigma[s])
      pi[s,g]     <- (2 * midpt[g] / (B * B)) * delta   # interval probability
      f[s,g]      <- p[s,g] * pi[s,g]
    }

    # Probability of capture (integrated across distance bands)
    pcap[s] <- sum(f[s,])

    # Normalize detection probabilities per band
    for (g in 1:nD) {
      fc[s,g] <- f[s,g] / pcap[s]
    }

    # Marginalized likelihood (fast mixing)
    ncap[s] ~ dpois(lambda[s] * pcap[s])

    # Abundance linear predictor
    log(lambda[s]) <-
        beta_int + beta_field * field[s] +
        beta_Y10 * Y10[s] + beta_Y12 * Y12[s] + beta_Y13 * Y13[s] +
        beta_Y19 * Y19[s] + beta_Y24 * Y24[s] + beta_Y25 * Y25[s] +
        beta_fieldY10 * field[s] * Y10[s] +
        beta_fieldY12 * field[s] * Y12[s] +
        beta_fieldY13 * field[s] * Y13[s] +
        beta_fieldY19 * field[s] * Y19[s] +
        beta_fieldY24 * field[s] * Y24[s] +
        beta_fieldY25 * field[s] * Y25[s]
  }

  #### ------------------------
  #### Derived parameters
  #### ------------------------

  # Mean per-survey abundance per site type × year (divided by 4 visits)
  lam.S.Y10 <- exp(beta_int + beta_field + beta_Y10 + beta_fieldY10) / 4
  lam.K.Y10 <- exp(beta_int            + beta_Y10) / 4

  lam.S.Y12 <- exp(beta_int + beta_field + beta_Y12 + beta_fieldY12) / 4
  lam.K.Y12 <- exp(beta_int            + beta_Y12) / 4

  lam.S.Y13 <- exp(beta_int + beta_field + beta_Y13 + beta_fieldY13) / 4
  lam.K.Y13 <- exp(beta_int            + beta_Y13) / 4

  lam.S.Y18 <- exp(beta_int + beta_field) / 4
  lam.K.Y18 <- exp(beta_int) / 4

  lam.S.Y19 <- exp(beta_int + beta_field + beta_Y19 + beta_fieldY19) / 4
  lam.K.Y19 <- exp(beta_int            + beta_Y19) / 4

  lam.S.Y24 <- exp(beta_int + beta_field + beta_Y24 + beta_fieldY24) / 4
  lam.K.Y24 <- exp(beta_int            + beta_Y24) / 4

  lam.S.Y25 <- exp(beta_int + beta_field + beta_Y25 + beta_fieldY25) / 4
  lam.K.Y25 <- exp(beta_int            + beta_Y25) / 4

  # Area of circular transect (ha)
  area <- 3.141592653589793 * B * B / 10000

  # Density (individuals per ha)
  D.S.Y10 <- lam.S.Y10 / area
  D.K.Y10 <- lam.K.Y10 / area

  D.S.Y12 <- lam.S.Y12 / area
  D.K.Y12 <- lam.K.Y12 / area

  D.S.Y13 <- lam.S.Y13 / area
  D.K.Y13 <- lam.K.Y13 / area

  D.S.Y18 <- lam.S.Y18 / area
  D.K.Y18 <- lam.K.Y18 / area

  D.S.Y19 <- lam.S.Y19 / area
  D.K.Y19 <- lam.K.Y19 / area

  D.S.Y24 <- lam.S.Y24 / area
  D.K.Y24 <- lam.K.Y24 / area

  D.S.Y25 <- lam.S.Y25 / area
  D.K.Y25 <- lam.K.Y25 / area
}
",fill=TRUE, file="JAGS/distance_JAGS_model_2025.txt")

# Inits
inits.fun()

# Params to save
params <- c("mu_alpha", "alpha_year_used", 
            "beta_int", "beta_field", "beta_Y10", "beta_Y12", 
            "beta_Y13", "beta_Y19", "beta_Y24",
            "beta_Y25", "beta_fieldY10",  "beta_fieldY12", 
            "beta_fieldY13", "beta_fieldY19",
            "beta_fieldY24", "beta_fieldY25",
            "D.S.Y10", "D.K.Y10", 
            "D.S.Y12", "D.K.Y12", 
            "D.S.Y13", "D.K.Y13",            
            "D.S.Y18", "D.K.Y18", 
            "D.S.Y19", "D.K.Y19", 
            "D.S.Y24", "D.K.Y24",
            "D.S.Y25", "D.K.Y25",
            "lam.S.Y10", "lam.K.Y10",
            "lam.S.Y12", "lam.K.Y12",
            "lam.S.Y13", "lam.K.Y13",
            "lam.S.Y18", "lam.K.Y18",
            "lam.S.Y19", "lam.K.Y19",
            "lam.S.Y24", "lam.K.Y24",
            "sigma")

# MCMC settings
ni <- 9000   ;   nb <- 100   ; na <- 3000;   nt <- 30   ;   nc <- 3

# Run JAGS and summarize posteriors
out <- jagsUI::jags(data = b, 
                    inits = inits.fun, 
                    parameters.to.save = params, 
                    model.file = "JAGS/distance_JAGS_model_2025.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
                    n.adapt = na,
                    n.cores = nc)
# sink("output/bobolink_2025.txt")
print(out, 2)
# sink()
# jagsUI:: traceplot(out)
# saveRDS(out, "output/bobo.psi.year.field.int.p.year.rds")




### PLOT


out <- readRDS("output/b.6years.int.model.p.year.rds")

# format data for plotting
b.plot <- data.frame(field = rep(c("S", "K"), 6),
                     year = c(2010, 2010, 2012, 2012, 2013, 2013, 
                              2018, 2018, 2019, 2019, 2024, 2024),
                     Nmed = c(median(out$sims.list$lam.S.Y10), median(out$sims.list$lam.K.Y10), 
                              median(out$sims.list$lam.S.Y12), median(out$sims.list$lam.K.Y12), 
                              median(out$sims.list$lam.S.Y13), median(out$sims.list$lam.K.Y13),
                              median(out$sims.list$lam.S.Y18), median(out$sims.list$lam.K.Y18), 
                              median(out$sims.list$lam.S.Y19), median(out$sims.list$lam.K.Y19), 
                              median(out$sims.list$lam.S.Y24), median(out$sims.list$lam.K.Y24)),
                     Nq2.5 = c(quantile(out$sims.list$lam.S.Y10, 0.025), 
                               quantile(out$sims.list$lam.K.Y10, 0.025), 
                               quantile(out$sims.list$lam.S.Y12, 0.025), 
                               quantile(out$sims.list$lam.K.Y12, 0.025), 
                               quantile(out$sims.list$lam.S.Y13, 0.025), 
                               quantile(out$sims.list$lam.K.Y13, 0.025),
                               quantile(out$sims.list$lam.S.Y18, 0.025), 
                               quantile(out$sims.list$lam.K.Y18, 0.025), 
                               quantile(out$sims.list$lam.S.Y19, 0.025), 
                               quantile(out$sims.list$lam.K.Y19, 0.025), 
                               quantile(out$sims.list$lam.S.Y24, 0.025), 
                               quantile(out$sims.list$lam.K.Y24, 0.025)),
                     Nq97.5 = c(quantile(out$sims.list$lam.S.Y10, 0.975), 
                                quantile(out$sims.list$lam.K.Y10, 0.975), 
                                quantile(out$sims.list$lam.S.Y12, 0.975), 
                                quantile(out$sims.list$lam.K.Y12, 0.975), 
                                quantile(out$sims.list$lam.S.Y13, 0.975), 
                                quantile(out$sims.list$lam.K.Y13, 0.975),
                                quantile(out$sims.list$lam.S.Y18, 0.975), 
                                quantile(out$sims.list$lam.K.Y18, 0.975), 
                                quantile(out$sims.list$lam.S.Y19, 0.975), 
                                quantile(out$sims.list$lam.K.Y19, 0.975), 
                                quantile(out$sims.list$lam.S.Y24, 0.975), 
                                quantile(out$sims.list$lam.K.Y24, 0.975)),           
                     Nq10 = c(quantile(out$sims.list$lam.S.Y10, 0.1), 
                              quantile(out$sims.list$lam.K.Y10, 0.1), 
                              quantile(out$sims.list$lam.S.Y12, 0.1), 
                              quantile(out$sims.list$lam.K.Y12, 0.1), 
                              quantile(out$sims.list$lam.S.Y13, 0.1), 
                              quantile(out$sims.list$lam.K.Y13, 0.1),
                              quantile(out$sims.list$lam.S.Y18, 0.1), 
                              quantile(out$sims.list$lam.K.Y18, 0.1), 
                              quantile(out$sims.list$lam.S.Y19, 0.1), 
                              quantile(out$sims.list$lam.K.Y19, 0.1), 
                              quantile(out$sims.list$lam.S.Y24, 0.1), 
                              quantile(out$sims.list$lam.K.Y24, 0.1)),
                     Nq90 = c(quantile(out$sims.list$lam.S.Y10, 0.9), 
                              quantile(out$sims.list$lam.K.Y10, 0.9), 
                              quantile(out$sims.list$lam.S.Y12, 0.9), 
                              quantile(out$sims.list$lam.K.Y12, 0.9), 
                              quantile(out$sims.list$lam.S.Y13, 0.9), 
                              quantile(out$sims.list$lam.K.Y13, 0.9),
                              quantile(out$sims.list$lam.S.Y18, 0.9), 
                              quantile(out$sims.list$lam.K.Y18, 0.9), 
                              quantile(out$sims.list$lam.S.Y19, 0.9), 
                              quantile(out$sims.list$lam.K.Y19, 0.9), 
                              quantile(out$sims.list$lam.S.Y24, 0.9), 
                              quantile(out$sims.list$lam.K.Y24, 0.9))) %>%
  mutate(Dmed = Nmed/pi,
         Dq2.5 = Nq2.5/pi,
         Dq97.5 = Nq97.5/pi,
         Dq10 = Nq10/pi,
         Dq90 = Nq90/pi)

(bobo <- ggplot(b.plot) +
    geom_point(aes(x = year, y = Dmed, color = field), size = 4,
               position = position_dodge(0.5)) +
    geom_errorbar(aes(x = year, ymin = Dq2.5, ymax = Dq97.5,
                      color = field),
                  width = 0,
                  position = position_dodge(0.5)) +
    geom_errorbar(aes(x = year, ymin = Dq10, ymax = Dq90,
                      color = field),
                  linewidth = 1,
                  width = 0,
                  position = position_dodge(0.5)) +
    geom_line(aes(x = year, y = Dmed, color = field),
              position = position_dodge(0.5)) +
    scale_color_manual(values = c("steelblue", "darkred")) +
    scale_x_continuous(breaks = seq(2010,2024, by = 4)) +
    ylim(0,6.1) +
    labs(y = "BOBO density (no./ha)", x = "", 
         color = "Field") +
    theme_bw() +
    theme(text = element_text(size = 13))
)

ggsave("figures/BOBO_6years.d.field.year.int.p.year.png", height = 4, width = 6, dpi = 400)



### PLOT RAW

# format raw count data
raw.plot <- b %>% filter(!is.na(dist)) %>% 
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
(raw.b <- ggplot(raw.plot) +
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

