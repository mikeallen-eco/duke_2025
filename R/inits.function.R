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
