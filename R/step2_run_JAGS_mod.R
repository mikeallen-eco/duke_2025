# function to run JAGS model for a species

library(jagsUI)

run_JAGS_mod <- function(jags_data, 
                         sp_name = "bobo",
                         mod,
                         ni = 130000, nb = 10000, 
                         na = 30000, nt = 1, nc = 3,
                         previous_model_run = NULL){
  
  # check if model output file exists first and only run if not
  if(!file.exists(file.path("output", paste0(sp_name, "_2025.rds")))){
    
  ### --- helper function to load inits
  Nst <- jags_data$ncap + 1
  nyears <- jags_data$nYears
  
    # If no previous model run provided → generic random initialization
    if (is.null(previous_model_run)) {message("No previous model run file detected. Using generic inits values...")}

    inits.fun <- function() {
        list(
          alpha_Int = rnorm(1, 4, 0.5),
          alpha_year_raw = rnorm(nyears, 0, 0.2),
          sd_year = runif(1, 0.1, 0.3),
          beta_Int = rnorm(1, 2, 0.5),
          beta_fieldS = rnorm(1, 0, 0.2),
          beta_Y10 = rnorm(1, 0, 0.2),
          beta_Y12 = rnorm(1, 0, 0.2),
          beta_Y13 = rnorm(1, 0, 0.2),
          beta_Y19 = rnorm(1, 0, 0.2),
          beta_Y24 = rnorm(1, 0, 0.2),
          beta_Y25 = rnorm(1, 0, 0.2),
          beta_fld10 = rnorm(1, 0, 0.2),
          beta_fld12 = rnorm(1, 0, 0.2),
          beta_fld13 = rnorm(1, 0, 0.2),
          beta_fld19 = rnorm(1, 0, 0.2),
          beta_fld24 = rnorm(1, 0, 0.2),
          beta_fld25 = rnorm(1, 0, 0.2),
          N = Nst
        )
    }
    
    inits_names <- names(inits.fun())
      
    # If previous model run provided → use posterior means as starting values
    if (!is.null(previous_model_run)) {
      message("Previous model run file detected. Using posterior means as inits...")
      inits.fun <- function() {
        mod <- readRDS(previous_model_run)
        
        # Check structure
        if (!"mean" %in% names(mod)) {
          stop("readRDS(previous_model_run) must have a 'mean' element containing posterior means.")
        }
        
        inits_list <- mod$mean
        
        # Drop derived parameters (e.g., those with 'lam', 'dev', or 'sig' in name)
        names_to_drop <- names(inits_list)[grepl("lam|dev|sig", names(inits_list))]
        inits_list <- inits_list[!(names(inits_list) %in% names_to_drop)]
        
        # Ensure N is properly initialized
        inits_list$N <- Nst
        
        # Return a plain list for JAGS
        as.list(inits_list)
      }
    }
  
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
              "sigma25",
              "pcap10", "pcap12", "pcap13",
              "pcap18", "pcap19", "pcap24",
              "pcap25")
  
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