# analyze grassland bird density at Duke Farms 2010-2025 with hierarchical distance sampling

source("R/setup.R")

### --- Step 1. Read in and format data

# read in data for each set of years
d10_13 <- format_2010_2013_data_for_analysis("data/Grassland_Birds_2010-2013_data.csv")
d18_19 <- format_2018_2019_data_for_analysis("data/Grassland_Birds_2018-2019_data.csv")
d24 <- format_2024_2025_data_for_analysis("data/Grassland_Birds_2024_data.csv")
d25 <- format_2024_2025_data_for_analysis("data/20250520_grassland_birds_2025_Allen/Grassland_Birds_2025_data.csv")
d <- bind_rows(d18_19,d24,d25)

# format JAGS data list for each species
b <- get_jags_data_for_distance_sampling(spname = "BOBO", dist_threshold = 100)

### --- Step 2. Edit & save JAGS model

# BUGS model specification for point transect data
cat("
model{
  # Priors
  alpha.int ~ dunif(-10,10)
  alpha ~ dunif(-5,5)
  alpha2 ~ dunif(-5,5)
  alpha3 ~ dunif(-5,5)
  alpha4 ~ dunif(-5,5)
  alpha5 ~ dunif(-5,5)
  beta0 ~ dunif(-5,5)
  beta1 ~ dunif(-5,5) # field covariate
  beta2 ~ dunif(-5,5) # Y10 covariate
  beta3 ~ dunif(-5,5) # Y12 covariate
  beta4 ~ dunif(-5,5) # Y13 covariate
  beta5 ~ dunif(-5,5) # Y19 covariate
  beta6 ~ dunif(-5,5) # Y24 covariate
  beta7 ~ dunif(-5,5) # interaction covariate
  beta8 ~ dunif(-5,5) # interaction covariate
  beta9 ~ dunif(-5,5) # interaction covariate
  beta10 ~ dunif(-5,5) # interaction covariate
  beta11 ~ dunif(-5,5) # interaction covariate

  for(i in 1:nind){
    dclass[i] ~ dcat(fc[site[i],]) # Part 1 of HM
  }
  for(s in 1:nsites){
    # Construct cell probabilities for nD distance bands
    for(g in 1:nD){                # midpt = mid-point of each band
      log(p[s,g]) <- -midpt[g] * midpt[g] / (2 * sigma[s] * sigma[s])
      pi[s,g] <- ((2 * midpt[g] ) / (B * B)) * delta # prob. per interval
      f[s,g] <- p[s,g] * pi[s,g]
      fc[s,g] <- f[s,g] / pcap[s]
    }
    pcap[s] <- sum(f[s,])           # Pr(capture): sum of rectangular areas

    ncap[s] ~ dbin(pcap[s], N[s])   # Part 2 of HM
    N[s] ~ dpois(lambda[s])         # Part 3 of HM
    log(lambda[s]) <- beta0 + beta1 * field[s] + 
        beta2 * Y10[s] + beta3 * Y12[s] +
        beta4 * Y13[s] + beta5 * Y19[s] + 
        beta6 * Y24[s] + 
        beta7 * field[s] * Y10[s] +
        beta8 * field[s] * Y12[s] +
        beta9 * field[s] * Y13[s] +
        beta10 * field[s] * Y19[s] +
        beta11 * field[s] * Y24[s] # linear model abundance
    log(sigma[s]) <- alpha0 + alpha1 * Y10[s] + alpha2 * Y12[s] + alpha3 * Y13[s] + 
                        alpha4 * Y19[s] +
                        alpha5 * Y24[s] # linear model detection
  }

  # Derived parameters
  lam.S.Y18 <- exp(beta0 + beta1)/4
  lam.K.Y18 <- exp(beta0)/4
  lam.S.Y10 <- exp(beta0 + beta1 + beta2 + beta7)/4
  lam.K.Y10 <- exp(beta0 + beta2)/4
  lam.S.Y12 <- exp(beta0 + beta1 + beta3 + beta8)/4
  lam.K.Y12 <- exp(beta0 + beta3)/4
  lam.S.Y13 <- exp(beta0 + beta1 + beta4 + beta9)/4
  lam.K.Y13 <- exp(beta0 + beta4)/4
  lam.S.Y19 <- exp(beta0 + beta1 + beta5 + beta10)/4
  lam.K.Y19 <- exp(beta0 + beta5)/4
  lam.S.Y24 <- exp(beta0 + beta1 + beta6 + beta11)/4
  lam.K.Y24 <- exp(beta0 + beta6)/4
  area.S <- 3.141*100*100/10000 # converted to ha [using 100 m as B was set to 103]
  area.K <- 3.141*100*100/10000 # converted to ha [using 100 m as B was set to 103]
  D.S.Y10 <- (lam.S.Y10)/area.S # mean per sampling period per ha
  D.K.Y10 <- (lam.K.Y10)/area.K # mean per sampling period per ha  
  D.S.Y12 <- (lam.S.Y12)/area.S # mean per sampling period per ha
  D.K.Y12 <- (lam.K.Y12)/area.K # mean per sampling period per ha  
  D.S.Y13 <- (lam.S.Y13)/area.S # mean per sampling period per ha
  D.K.Y13 <- (lam.K.Y13)/area.K # mean per sampling period per ha  
  D.S.Y18 <- (lam.S.Y18)/area.S # mean per sampling period per ha
  D.K.Y18 <- (lam.K.Y18)/area.K # mean per sampling period per ha  
  D.S.Y19 <- (lam.S.Y19)/area.S # mean per sampling period per ha
  D.K.Y19 <- (lam.K.Y19)/area.K # mean per sampling period per ha  
  D.S.Y24 <- (lam.S.Y24)/area.S # mean per sampling period per ha
  D.K.Y24 <- (lam.K.Y24)/area.K # mean per sampling period per ha
  sigma10 <- sigma[1]
  sigma12 <- sigma[20]
  sigma13 <- sigma[38]
  sigma18 <- sigma[57]
  sigma19 <- sigma[76]
  sigma24 <- sigma[95]
}
",fill=TRUE, file="scripts/distance.N.field.6year.int.p.6year.txt")


# Inits
Nst <- ncap + 1
inits.fun <- function(){list(alpha0=3.8, alpha1=0, alpha2=0,
                             beta0=1.8, beta1=0.5, beta2 = 1, beta3 = .8, 
                             beta4 = 0, beta5 = -1, beta6 = 0,
                             beta7 = -.6, beta8 = 0, beta9 = 0.5,
                             beta10 = 0, beta11 = 0, N=Nst)}

# Params to save
params <- c("alpha0", "alpha1", "alpha2",
            "alpha3", "alpha4", "alpha5",
            "beta0", "beta1",
            "beta2", "beta3", 
            "beta4", "beta5",
            "beta6", "beta7",
            "beta8", "beta9",
            "beta10", "beta11",
            "D.S.Y10", "D.K.Y10", 
            "D.S.Y12", "D.K.Y12", 
            "D.S.Y13", "D.K.Y13",            
            "D.S.Y18", "D.K.Y18", 
            "D.S.Y19", "D.K.Y19", 
            "D.S.Y24", "D.K.Y24",
            "lam.S.Y10", "lam.K.Y10",
            "lam.S.Y12", "lam.K.Y12",
            "lam.S.Y13", "lam.K.Y13",
            "lam.S.Y18", "lam.K.Y18",
            "lam.S.Y19", "lam.K.Y19",
            "lam.S.Y24", "lam.K.Y24",
            "sigma10", "sigma12", "sigma13",
            "sigma18", "sigma19", "sigma24") 

# MCMC settings
ni <- 90000   ;   nb <- 0   ; na <- 30000;   nt <- 30   ;   nc <- 3

# Run BUGS (~ 2.3 min) and summarize posteriors
out <- jagsUI::jags(data = b.jags.data, 
                    inits = inits.fun, 
                    parameters.to.save = params, 
                    model.file = "scripts/distance.N.field.6year.int.p.6year.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
                    n.adapt = na,
                    n.cores = 3)
# sink("output/b.6years.model.int.p.year.txt")
print(out, 2)
# sink()
# jagsUI:: traceplot(out)
# saveRDS(out, "output/b.6years.int.model.p.year.rds")




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

