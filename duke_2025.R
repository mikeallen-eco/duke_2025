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

### --- Step 2. Edit & save JAGS model

# JAGS model specification for point transect data
cat("
model{
  # Priors
  alpha_Int ~ dunif(-5,5)
  alpha_Y10 ~ dunif(-5,5)
  alpha_Y12 ~ dunif(-5,5)
  alpha_Y13 ~ dunif(-5,5)
  alpha_Y19 ~ dunif(-5,5)
  alpha_Y24 ~ dunif(-5,5)
  alpha_Y25 ~ dunif(-5,5)
  beta_Int ~ dunif(-5,5)
  beta_fieldS ~ dunif(-5,5) # field covariate
  beta_Y10 ~ dunif(-5,5) # Y10 covariate
  beta_Y12 ~ dunif(-5,5) # Y12 covariate
  beta_Y13 ~ dunif(-5,5) # Y13 covariate
  beta_Y19 ~ dunif(-5,5) # Y19 covariate
  beta_Y24 ~ dunif(-5,5) # Y24 covariate
  beta_Y25 ~ dunif(-5,5) # Y25 covariate
  beta_fld10 ~ dunif(-5,5) # interaction covariate
  beta_fld12 ~ dunif(-5,5) # interaction covariate
  beta_fld13 ~ dunif(-5,5) # interaction covariate
  beta_fld19 ~ dunif(-5,5) # interaction covariate
  beta_fld24 ~ dunif(-5,5) # interaction covariate
  beta_fld25 ~ dunif(-5,5) # interaction covariate

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
    
    # linear model of abundance
    log(lambda[s]) <- beta_Int + beta_fieldS * field[s] + 
        beta_Y10 * Y10[s] + beta_Y12 * Y12[s] +
        beta_Y13 * Y13[s] + beta_Y19 * Y19[s] + 
        beta_Y24 * Y24[s] + beta_Y25 * Y25[s] +
        beta_fld10 * field[s] * Y10[s] +
        beta_fld12 * field[s] * Y12[s] +
        beta_fld13 * field[s] * Y13[s] +
        beta_fld19 * field[s] * Y19[s] +
        beta_fld24 * field[s] * Y24[s] +
        beta_fld25 * field[s] * Y25[s]
    
    # linear model of detection
    log(sigma[s]) <- alpha_Int + alpha_Y10 * Y10[s] + alpha_Y12 * Y12[s] +
                      alpha_Y13 * Y13[s] + alpha_Y19 * Y19[s] +
                        alpha_Y24 * Y24[s] + alpha_Y25 * Y25[s]
  }

  # Derived parameters
  
  # lambda: mean count per sample
  lam.S.Y10 <- exp(beta_Int + beta_fieldS + beta_Y10 + beta_fld10)/4
  lam.K.Y10 <- exp(beta_Int + beta_Y10)/4
  lam.S.Y12 <- exp(beta_Int + beta_fieldS + beta_Y12 + beta_fld12)/4
  lam.K.Y12 <- exp(beta_Int + beta_Y12)/4
  lam.S.Y13 <- exp(beta_Int + beta_fieldS + beta_Y13 + beta_fld13)/4
  lam.K.Y13 <- exp(beta_Int + beta_Y13)/4
  lam.S.Y18 <- exp(beta_Int + beta_fieldS)/4
  lam.K.Y18 <- exp(beta_Int)/4
  lam.S.Y19 <- exp(beta_Int + beta_fieldS + beta_Y19 + beta_fld19)/4
  lam.K.Y19 <- exp(beta_Int + beta_Y19)/4
  lam.S.Y24 <- exp(beta_Int + beta_fieldS + beta_Y24 + beta_fld24)/4
  lam.K.Y24 <- exp(beta_Int + beta_Y24)/4
  lam.S.Y25 <- exp(beta_Int + beta_fieldS + beta_Y25 + beta_fld25)/4
  lam.K.Y25 <- exp(beta_Int + beta_Y25)/4
  area.S <- 3.141*100*100/10000 # converted to ha
  area.K <- 3.141*100*100/10000 # converted to ha
  
  # density: mean count per sampling period per ha
  D.S.Y10 <- (lam.S.Y10)/area.S
  D.K.Y10 <- (lam.K.Y10)/area.K  
  D.S.Y12 <- (lam.S.Y12)/area.S
  D.K.Y12 <- (lam.K.Y12)/area.K  
  D.S.Y13 <- (lam.S.Y13)/area.S
  D.K.Y13 <- (lam.K.Y13)/area.K  
  D.S.Y18 <- (lam.S.Y18)/area.S
  D.K.Y18 <- (lam.K.Y18)/area.K 
  D.S.Y19 <- (lam.S.Y19)/area.S
  D.K.Y19 <- (lam.K.Y19)/area.K  
  D.S.Y24 <- (lam.S.Y24)/area.S
  D.K.Y24 <- (lam.K.Y24)/area.K
  D.S.Y25 <- (lam.S.Y25)/area.S
  D.K.Y25 <- (lam.K.Y25)/area.K
  
  # annual sigma for half-normal detection function
  sigma10 <- sigma[1]
  sigma12 <- sigma[20]
  sigma13 <- sigma[38]
  sigma18 <- sigma[57]
  sigma19 <- sigma[76]
  sigma24 <- sigma[95]
  sigma25 <- sigma[114]
}
",fill=TRUE, file="JAGS/distance_model_2025.txt")

# Params to save
params <- c("alpha_Int", "alpha_Y10", "alpha_Y12",
            "alpha_Y13", "alpha_Y19", "alpha_Y24",
            "alpha_Y25",
            "beta_Int", "beta_fieldS",
            "beta_Y10", "beta_Y12", 
            "beta_Y13", "beta_Y19",
            "beta_Y24", "beta_Y25",
            "beta_fld10", "beta_fld12",
            "beta_fld13", "beta_fld19",
            "beta_fld24", "beta_fld25",
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
            "lam.S.Y25", "lam.K.Y25",
            "sigma10", "sigma12", "sigma13",
            "sigma18", "sigma19", "sigma24",
            "sigma25")

# MCMC settings
ni <- 130000   ;   nb <- 10000; na <- 30000;   nt <- 30   ;   nc <- 3

# --- BOBOLINK

# Inits
bNst <- b$ncap + 1
inits.fun <- function(){list(alpha_Int=4, alpha_Y10=0, alpha_Y12=0,
                             alpha_Y13=0.5, alpha_Y19=0, alpha_Y24=0, 
                             alpha_Y25=0,
                             beta_Int=2, beta_fieldS = 0, beta_Y10 = 0, 
                             beta_Y12 = 0, beta_Y13 = 0, beta_Y19 = 0, 
                             beta_Y24 = 0.5, beta_Y25 = 0.5, beta_fld10 = 0, 
                             beta_fld12 = 0, beta_fld13 = 0, beta_fld19 = 0, 
                             beta_fld24 = 0, beta_fld25 = 0, N=bNst)}

# Run JAGS and summarize posteriors
bobo_out <- jagsUI::jags(data = b, 
                    inits = inits.fun, 
                    parameters.to.save = params, 
                    model.file = "JAGS/distance_model_2025.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
                    n.adapt = na,
                    n.cores = nc)
sink("output/bobo_2025.txt")
print(bobo_out, 2)
sink()
# jagsUI:: traceplot(bobo_out)
saveRDS(bobo_out, "output/bobo_2025.rds")

# --- GRASSHOPPER SPARROW

# Inits
gNst <- g$ncap + 1
inits.fun <- function(){list(alpha_Int=4, alpha_Y10=0, alpha_Y12=0,
                             alpha_Y13=0.5, alpha_Y19=0, alpha_Y24=0, 
                             alpha_Y25=0,
                             beta_Int=2, beta_fieldS = 0, beta_Y10 = 0, 
                             beta_Y12 = 0, beta_Y13 = 0, beta_Y19 = 0, 
                             beta_Y24 = 0.5, beta_Y25 = 0.5, beta_fld10 = 0, 
                             beta_fld12 = 0, beta_fld13 = 0, beta_fld19 = 0, 
                             beta_fld24 = 0, beta_fld25 = 0, N=gNst)}

# Run JAGS and summarize posteriors
grsp_out <- jagsUI::jags(data = g, 
                         inits = inits.fun, 
                         parameters.to.save = params, 
                         model.file = "JAGS/distance_model_2025.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
                         n.adapt = na,
                         n.cores = nc)
sink("output/grsp_2025.txt")
print(grsp_out, 2)
sink()
# jagsUI:: traceplot(grsp_out)
saveRDS(grsp_out, "output/grsp_2025.rds")

# --- EASTERN MEADOWLARK

# Inits
eNst <- e$ncap + 1
inits.fun <- function(){list(alpha_Int=4, alpha_Y10=0, alpha_Y12=0,
                             alpha_Y13=0.5, alpha_Y19=0, alpha_Y24=0, 
                             alpha_Y25=0,
                             beta_Int=2, beta_fieldS = 0, beta_Y10 = 0, 
                             beta_Y12 = 0, beta_Y13 = 0, beta_Y19 = 0, 
                             beta_Y24 = 0.5, beta_Y25 = 0.5, beta_fld10 = 0, 
                             beta_fld12 = 0, beta_fld13 = 0, beta_fld19 = 0, 
                             beta_fld24 = 0, beta_fld25 = 0, N=eNst)}

# Run JAGS and summarize posteriors
eame_out <- jagsUI::jags(data = e, 
                         inits = inits.fun, 
                         parameters.to.save = params, 
                         model.file = "JAGS/distance_model_2025.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
                         n.adapt = na,
                         n.cores = nc)
sink("output/eame_2025.txt")
print(eame_out, 2)
sink()
# jagsUI:: traceplot(eame_out)
saveRDS(eame_out, "output/eame_2025.rds")

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

