# function to run JAGS model for a species

library(jagsUI)

run_JAGS_mod <- function(jags_data, 
                         sp_name = "bobo",
                         mod = "JAGS/distance_model_2025.txt",
                         ni = 130000, nb = 10000, 
                         na = 30000, nt = 30, nc = 3){
  
  # check if model output file exists first and only run if not
  if(!file.exists(file.path("output", paste0(sp_name, "_2025.rds")))){
    
  ### --- helper function to load inits
  
  inits.fun <- function() {
    list(
      # --- Abundance model ---
      beta_int       = rnorm(1, 0, 1),     # mean log-abundance
      beta_field     = rnorm(1, 0, 0.5),   # modest site effect
      beta_Y10       = rnorm(1, 0, 0.5),
      beta_Y12       = rnorm(1, 0, 0.5),
      beta_Y13       = rnorm(1, 0, 0.5),
      beta_Y19       = rnorm(1, 0, 0.5),
      beta_Y24       = rnorm(1, 0, 0.5),
      beta_Y25       = rnorm(1, 0, 0.5),
      
      beta_fieldY10  = rnorm(1, 0, 0.5),
      beta_fieldY12  = rnorm(1, 0, 0.5),
      beta_fieldY13  = rnorm(1, 0, 0.5),
      beta_fieldY19  = rnorm(1, 0, 0.5),
      beta_fieldY24  = rnorm(1, 0, 0.5),
      beta_fieldY25  = rnorm(1, 0, 0.5),
      
      # --- Detection model ---
      mu_alpha       = rnorm(1, 3.5, 0.2),      # mean log(sigma) ≈ 3.5
      sigma_alpha    = runif(1, 0.5, 1.5),      # detection SD
      alpha_year     = rnorm(nyears, 3.5, 0.5)  # year-level random effects
    )
  }
  
  ### --- define model parameters to track
  
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
  
  ### --- define initial parameter states
  
  Nst <- jags_data$ncap + 1
  inits.fun <- function(){list(alpha_Int=4, alpha_Y10=0, alpha_Y12=0,
                               alpha_Y13=0.5, alpha_Y19=0, alpha_Y24=0, 
                               alpha_Y25=0,
                               beta_Int=2, beta_fieldS = 0, beta_Y10 = 0, 
                               beta_Y12 = 0, beta_Y13 = 0, beta_Y19 = 0, 
                               beta_Y24 = 0.5, beta_Y25 = 0.5, beta_fld10 = 0, 
                               beta_fld12 = 0, beta_fld13 = 0, beta_fld19 = 0, 
                               beta_fld24 = 0, beta_fld25 = 0, N=bNst)}
  
  ### --- Run JAGS and save output
  
  out <- jagsUI::jags(data = jags_data, 
                      inits = inits.fun, 
                      parameters.to.save = params, 
                      model.file = "JAGS/distance_model_2025.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
                      n.adapt = na,
                      n.cores = nc)
  sink(file.path("output", paste0(sp_name, "_2025.txt")))
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