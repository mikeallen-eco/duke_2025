
get_alpha_priors <- function(mod = "output/bobo_2025.rds"){

# read in JAGS model object
out <- readRDS(mod)

# combine posteriors for years with alpha slope estimates
sims <- c(out$sims.list$alpha_Y19, out$sims.list$alpha_Y24, out$sims.list$alpha_Y25)

# get posterior mean and SD for those pooled posteriors
mu <- mean(sims)
s <- sd(sims)
tau <- 1 / (s^2)

return(list(alpha_mu = mu,
            alpha_sd = s,
            alpha_tau = tau))

}
