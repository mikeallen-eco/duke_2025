get_HDS_density <- function(mod){
  
  model <- readRDS(mod)
  dens_S_25 <- model$sims.list$lam.S.Y25 / pi
  dens_K_25 <- model$sims.list$lam.K.Y25 / pi
  dens_S_24 <- model$sims.list$lam.S.Y24 / pi
  dens_K_24 <- model$sims.list$lam.K.Y24 / pi
  
  chg_S_25_24 <- 100*((dens_S_25 - dens_S_24) / (dens_S_24))
  chg_K_25_24 <- 100*((dens_K_25 - dens_K_24) / (dens_K_24))
  
  list(S_med = median(dens_S_25),
       S_q2.5 = quantile(dens_S_25, 0.025),
       S_q97.5 = quantile(dens_S_25, 0.975),
       K_med = median(dens_K_25),
       K_q2.5 = quantile(dens_K_25, 0.025),
       K_q97.5 = quantile(dens_K_25, 0.975),
       chg_S_med = median(chg_S_25_24),
       chg_S_q2.5 = quantile(chg_S_25_24, 0.025),
       chg_S_q97.5 = quantile(chg_S_25_24, 0.975),
       chg_K_med = median(chg_K_25_24),
       chg_K_q2.5 = quantile(chg_K_25_24, 0.025),
       chg_K_q97.5 = quantile(chg_K_25_24, 0.975))
  
  
}