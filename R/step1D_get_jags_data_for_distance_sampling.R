# format JAGS data list for a species for distance sampling

library(tidyr)

get_jags_data_for_distance_sampling <- function(df = d, spname = "BOBO", dist_threshold = 100){
  
  # get info for each point-year (i.e., a 'survey' in this analysis)
  survey_info <- df %>%
    select(pt, field, year) %>%
    distinct() %>%
    arrange(pt, year)
  
  # extract all observations of the species (years 2018+ when distance data is available)
  sp_data_later_years <- df %>%
    filter(!year %in% c("Y10", "Y12", "Y13")) %>%
    filter(species %in% spname, dist < (dist_threshold + 1))

  # extract all observations of the species (years 2010-2013 when distance data not available)
  sp_data_early_years <- df %>%
    filter(year %in% c("Y10", "Y12", "Y13")) %>%
    filter(species %in% spname)
    # filter(species %in% spname, dist < (dist_threshold + 1))
  
  # impute distances for 2010-2013 (same bin probabilities as 2018+)
  bins <- cut(sp_data_later_years$dist, breaks = seq(0, 100, by = 5), include.lowest = TRUE, right = FALSE)
  bin_probs <- prop.table(table(bins)) # Compute bin probabilities
  n_missing <- length(sp_data_early_years$dist)
  set.seed(102) # Sample according to observed proportions
  sampled_bins <- sample(names(bin_probs), size = n_missing, replace = TRUE, prob = bin_probs)
  bin_edges <- as.data.frame(do.call(rbind, strsplit(gsub("\\[|\\)|\\]", "", sampled_bins), ",")))
  colnames(bin_edges) <- c("lower", "upper")
  bin_edges <- bin_edges |> dplyr::mutate(lower = as.numeric(lower), upper = as.numeric(upper))
  sp_data_early_years$dist <- runif(n_missing, min = bin_edges$lower, max = bin_edges$upper)
  sp_data_early_years$dist <- 9999
  
  # recombine 2010-2013 data to 2018+ data
  sp_data <- sp_data_early_years %>%
    bind_rows(sp_data_later_years) %>%
    select(-distcat)
  
  # get names of surveys that had zero observations of the species
  zero_surveys <- setdiff(unique(survey_info$pt), unique(sp_data$pt))
  
  # add zero surveys back in and format covariates
  sp_data_w_zeros <- sp_data %>%
    bind_rows(survey_info %>% filter(pt %in% zero_surveys)) %>%
    select(field, year, pt, dist, num) %>%
    arrange(pt) %>%
    mutate(
      ptn = as.numeric(as.factor(pt)),
      dist.na = ifelse(is.na(dist), NA, 1),
      fld = as.numeric(as.factor(field)) - 1
    ) %>%
    select(field, year, pt, dist, ptn, dist.na, fld, num)
  
  # expand df so that each observation has one row
  sp_data_expanded <- sp_data_w_zeros %>%
    mutate(num_orig = num, num = ifelse(is.na(num), 0, num)) %>%
    mutate(num = ifelse(num == 0, 1, num)) %>%
    tidyr::uncount(weights = num, .remove = FALSE) %>%
    mutate(num = ifelse(is.na(num_orig), 0, 1)) %>%
    select(-num_orig)
  
  # make field covariate (1 row for each site; 1 = skeet, 0 = kaufman)
  site.covs <- sp_data_expanded %>%
    select(ptn, year, field) %>%
    distinct() %>%
    mutate(
      fld = as.numeric(as.factor(field)) - 1,
      Y10 = as.numeric(ifelse(year %in% "Y10", 1, 0)),
      Y12 = as.numeric(ifelse(year %in% "Y12", 1, 0)),
      Y13 = as.numeric(ifelse(year %in% "Y13", 1, 0)),
      Y18 = as.numeric(ifelse(year %in% "Y18", 1, 0)),
      Y19 = as.numeric(ifelse(year %in% "Y19", 1, 0)),
      Y24 = as.numeric(ifelse(year %in% "Y24", 1, 0)),
      Y25 = as.numeric(ifelse(year %in% "Y25", 1, 0))
    )
  
  # Prepare final data for JAGS
  ncap <- table(sp_data_expanded$ptn)            # ncap = 1 if no individuals captured
  sites0 <- sp_data_expanded[is.na(sp_data_expanded[, "dist.na"]), ]$ptn # sites where nothing was seen
  ncap[as.character(sites0)] <- 0    # Fill in 0 for sites with no detections
  ncap <- as.vector(ncap)            # Number of individuals detected per site
  B <- 103 # max radius distance
  
  # Other data
  site <- sp_data_expanded[!is.na(sp_data_expanded$dist.na),]$ptn   # Site ID of each observation
  delta <-  5                       # Distance bin width for rect. approx.
  midpt <- seq(delta/2, B, delta)    # Make mid-points and chop up data
  dclass <- sp_data_expanded$dist %/% delta + 1   # Convert distance to distance category
  nD <- length(midpt)               # Number of distance intervals
  dclass <- dclass[!is.na(sp_data_expanded[,"dist.na"])] # Observed categorical observations
  # special line to recode years 2010-2013 as dclass = NA (= unobserved and informed solely by the prior)
  dclass[dclass==2000] <- NA
  nind <- length(dclass)             # Total number of individuals detected
  nsites <- length(unique(site)) + length(unique(sites0))
  
  # define species number code for species-specific Year 2010-2013 priors
  if(spname %in% "BOBO"){species <- 1}
  if(spname %in% "GRSP"){species <- 2}
  if(spname %in% "EAME"){species <- 3}
  
  # Bundle and summarize data set
  jags.data <- list(nsites=nsites, nind=nind, B=B, nD=nD, midpt=midpt,
                    delta=delta, ncap=ncap, field=site.covs$fld, 
                    Y10 = site.covs$Y10, Y12 = site.covs$Y12, Y13 = site.covs$Y13,
                    Y19 = site.covs$Y19, Y24 = site.covs$Y24, Y25 = site.covs$Y25,
                    dclass=dclass, site=site, species=species)
  
  return(jags.data)
  
}