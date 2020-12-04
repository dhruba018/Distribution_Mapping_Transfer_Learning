## Matching function for mapping one distribution into another i.e., 
## distribution matching 
## 
## Dependency: stats 
## Dependency_own: lambda_functions 
################################################################################

match_func <- function(knots, vals, new_vals, lims) {
    
    source("lambda_functions.R")
    
    ## Limits for function inputs...
    if (missing(lims)) 
        lims <- c(min(knots), max(knots))
    
    ## Inverse CDF mapping...
    map <- stats::approxfun(x = vals, y = knots, yleft = lims[1], yright = lims[2], method = "linear", ties = "ordered", rule = 2)
    
    ## Get matched values...
    new_knots <- confined(map(new_vals), lims)
    new_knots

}