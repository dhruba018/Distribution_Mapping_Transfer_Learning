dist.match.trans.learn <- function(target.set, source.set, method = "hist", size = 1e3, seed = NULL, pred.opt = FALSE) {
    
    ## Initial check...
    if (ncol(target.set[["X"]]) != ncol(source.set[["X"]]))
        stop("Source and target set covariates must have the same number of features!")
    
    if (nrow(source.set[["X"]]) != length(source.set[["y"]]))
        stop("Source sets must have the same number of samples!")
    
    out01 <- function(y) (min(y) < 0) & (max(y) > 1)
    if (out01(target.set[["y"]]) | out01(source.set[["y"]]))
        stop("Response data must be normalized in [0, 1].")
    
    
    ## Load function files...
    source("lambda_functions.R")
    source("get_dist_est.R")
    source("dist_match.R")
    source("RF_predict.R")
    
    
    ######## MAIN ##############################################################
    
    ## Define datasets...
    X1 <- norm.data(target.set[["X"]]);      y1 <- target.set[["y"]]
    X2 <- norm.data(source.set[["X"]]);      y2 <- source.set[["y"]]
    n.feat <- ncol(X1);                      data.lims <- c(0, 1)
    
    
    ## Distribution matching for predictors...
    X2.map <- lapply(1:n.feat, function(j) { 
        dist_match(X1[, j], ref = X2[, j], match_method = method, samp_size = size, lims = data.lims)
        })
    X2.map <- as.data.frame(X2.map);    dimnames(X2.map) <- dimnames(X1)
    # rownames(X2.map) <- rownames(X1);   colnames(X2.map) <- colnames(X1)
    
    
    ## Perform prediction & map back to original space...
    y2.pred.map <- RF_predict(x_train = X2, y_train = y2, x_test = X2.map, y_lims = data.lims, 
                              n_tree = 200, m_try = 0.4, random_seed = seed)
    y2.dist <- get_dist_est(y2, sample_size = size, x_range = "unit", dist_method = method, grid_size = 1e3)
    
    y1.pred <- dist_match(y2.pred.map, ref = y1, src_cdf = y2.dist, match_method = method, samp_size = size, lims = data.lims)
    y1.pred <- confined(y1.pred, lims = data.lims);      names(y1.pred) <- names(y2.pred.map)
    
    
    ## Return output objects...
    if (pred.opt) {
        return( list("mapped" = y1.pred, "unmapped" = y2.pred.map) )
    } else {
        return( y1.pred )
    }

}