## Function for performing distribution matching between source and reference 
## datasets- extended from the idea of histogram matching (discrete) to matching
## continuous density estimates 
## 
## Dependency: stats, ks 
## Dependency_own: lambda_functions 
################################################################################

dist_match <- function(src, ref, src_cdf, ref_cdf, lims, match_method = "hist", samp_size = 1e6) {
  
  source("get_dist_est.R")
  source("match_func.R")
  
  ## Get distributions...
  if (missing(ref_cdf)) 
    ref_cdf <- get_dist_est(ref, sample_size = samp_size, x_range = "unit", dist_method = match_method, grid_size = 1e3)
  
  if (missing(src_cdf))  
    src_cdf <- get_dist_est(src, sample_size = samp.size, x_range = "unit", dist_method = match_method, grid_size = 1e3)
  
  
  ## Mapping parameters...
  match_method <- tolower(match_method)
  if (match_method == "hist") {                           # Using histogram
    kn_vals <- knots(ref_cdf);              fn_vals <- ref_cdf(kn_vals)
    vals_to_match <- src_cdf(src)
  } 
  
  else if (match_method == "dens") {                      # Using kernel density
    kn_vals <- ref_cdf$eval.points;         fn_vals <- ref_cdf$estimate
    vals_to_match <- predict(src_cdf, x = src)
  } 
  
  else {
    stop("Invalid option for estimating distribution! Please use 'hist' for histogram or 'dens' for kernel density.")
  }
  
  ## Perform mapping...
  if (missing(lims))
    lims <- c(min(src), max(src))
  
  matched <- match_func(knots = kn_vals, vals = fn_vals, new_vals = vals_to_match, lims)
  matched
  
}