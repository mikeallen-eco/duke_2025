# function to run JAGS model for a species

library(jagsUI)

run_JAGS_mod <- function(jags_data, 
                         sp_name = "bobo",
                         mod,
                         ni = 130000, nb = 10000, 
                         na = 30000, nt = 1, nc = 3){
  
  # check if model output file exists first and only run if not
  if(!file.exists(file.path("output", paste0(sp_name, "_2025.rds")))){
    
  ### --- helper function to load inits
  Nst <- jags_data$ncap + 1
  nyears <- jags_data$nYears
  # inits.fun <- function() {
  #   list(
  #     # --- Abundance model ---
  #     beta_Int       = rnorm(1, 2, .1),     # mean log-abundance
  #     beta_year      = rnorm(nyears, 0, 0.1),
  #     beta_fieldS     = rnorm(nyears, 1, 0.5),   # modest site effect
  #     
  #     # --- Detection model ---
  #     alpha_Int       = rnorm(1, 4, 0.1),      # mean log(sigma) ≈ 3.5
  #     sd_year         = 0.15,      # detection SD
  #     alpha_year      = rnorm(nyears, 0, 0.1),  # year-level random effects
  #   
  #     N=Nst)
  # }
  
  ### --- define initial parameter states
  
  inits.fun <- function(){list(alpha_Int=4, alpha_Y10=0, alpha_Y12=0,
                               alpha_Y13=0.5, alpha_Y19=0, alpha_Y24=0,
                               alpha_Y25=0, alpha_year=rep(0,nyears), sd_year = 0.2,
                               beta_Int=2, beta_fieldS = 0, beta_Y10 = 0,
                               beta_Y12 = 0, beta_Y13 = 0, beta_Y19 = 0,
                               beta_Y24 = 0, beta_Y25 = 0, beta_fld10 = 0,
                               beta_fld12 = 0, beta_fld13 = 0, beta_fld19 = 0,
                               beta_fld24 = 0, beta_fld25 = 0, N=Nst)}
  
  ### --- define model parameters to track
  
  params <- c("alpha_Int", "alpha_Y10", "alpha_Y12",
              "alpha_Y13", "alpha_Y19", "alpha_Y24",
              "alpha_Y25", "alpha_Y10_Y13", "alpha_year",
              "sd_year",
              "beta_Int", "beta_year", "beta_fieldS",
              "beta_Y10", "beta_Y12", 
              "beta_Y13", "beta_Y19",
              "beta_Y24", "beta_Y25",
              "beta_fld10", "beta_fld12",
              "beta_fld13", "beta_fld19",
              "beta_fld24", "beta_fld25",
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
  
  ### --- Run JAGS and save output
  
  out <- jagsUI::jags(data = jags_data, 
                      inits = inits.fun, 
                      parameters.to.save = params, 
                      model.file = mod,
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
                      n.adapt = na,
                      n.cores = nc)
  sink(file.path("output", paste0(sp_name, "_2025.txt")))
  print(out)
  sink()
  
  saveRDS(out, file.path("output", paste0(sp_name, "_2025.rds")))
  
  message("Done. Model output for ", sp_name, " saved in ", 
          file.path("output", paste0(sp_name, "_2025.rds")))
}else{
  message("Model file already exists. Returning saved model: ",
          file.path("output", paste0(sp_name, "_2025.rds")))
  out <- readRDS(file.path("output", paste0(sp_name, "_2025.rds")))
  }
  return(out)
  
}