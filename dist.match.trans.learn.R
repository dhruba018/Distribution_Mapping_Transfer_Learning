dist.match.trans.learn <- function(target.set, source.set, method = "hist", size = 1e3, seed = NULL, pred.opt = FALSE){
    
    ## Initial check...
    if (ncol(target.set[["X"]]) != ncol(source.set[["X"]]))
        stop("Source and target set covariates must have the same number of features!")
    
    if (nrow(source.set[["X"]]) != length(source.set[["y"]]))
        stop("Source sets must have the same number of samples!")
    
    out01 <- function(y) (min(y) < 0) & (max(y) > 1)
    if (out01(target.set[["y"]]) | out01(source.set[["y"]]))
        stop("Response data must be normalized in [0, 1].")
    
    
    ## Anonymous functions...
    norm01    <- function(z) { z <- if (min(z) > 0) z - min(z);   z <- z / max(z);  z }
    norm.data <- function(df) as.data.frame(apply(df, MARGIN = 2, norm01))
    # zscore    <- function(df) as.data.frame(apply(df, MARGIN = 2, scale))
    conf.lims <- function(y, lims = c(0, 1)){ y[y < lims[1]] <- lims[1];    y[y > lims[2]] <- lims[2];   y }
    
    
    ## Estimate distributions...
    get.cum.dist <- function(x, sample.size = 1e6, size.tol = 1e3, dist.method = "hist", grid.size = 1e4, ...) {
        dist.method <- tolower(dist.method)
        
        ## Bootstrapping...
        xx <- if (length(x) < size.tol) sample(x, size = sample.size, replace = TRUE) else x
        
        ## Cumulative distribution...
        if (dist.method == "hist") {                # Use histogram
            x.cdf <- ecdf(xx, ...)
        } 
        
        else if (dist.method == "dens") {           # Use kernel density
            if (!require(ks)) { library(ks) }       # Load package
            
            bw <- hscv(x, nstage = 2, binned = TRUE, bgridsize = grid.size * 10)
            x.cdf <- kcde(xx, h = bw, binned = TRUE, bgridsize = grid.size, xmin = 0, xmax = 1, ...)
        } 
        
        else {
            stop("Invalid option for estimating distribution! Please use 'hist' for histogram or 'dens' for kernel density.")
        }
        
        x.cdf
    }
    
    
    ## Matching function...
    match.func <- function(func.knots, func.vals, new.knots, knot.lims) {
        if (missing(knot.lims)) 
            knot.lims <- c(min(func.knots), max(func.knots))
        
        ## Inverse cdf mapping...
        map.func <- approxfun(x = func.vals, y = func.knots, yleft = knot.lims[1], yright = knot.lims[2], 
                              method = "linear", ties = "ordered", rule = 2)
        
        ## Get matched values...
        map.vals <- conf.lims(map.func(new.knots), lims = knot.lims)
        map.vals
    }
    
    
    ## Transfer learning function...
    dist.match <- function(src, ref, src.cdf, ref.cdf, lims, match.method = "hist", samp.size = 1e6) {
        ## Get distributions...
        if (missing(ref.cdf))
            ref.cdf <- get.cum.dist(ref, sample.size = samp.size, dist.method = match.method, grid.size = 1e3)
        
        if (missing(src.cdf)) 
            src.cdf <- get.cum.dist(src, sample.size = samp.size, dist.method = match.method, grid.size = 1e3)
        
        
        ## Mapping parameters...
        match.method <- tolower(match.method)
        if (match.method == "hist") {                           # Using histogram
            kn.vals <- knots(ref.cdf);              fn.vals <- ref.cdf(kn.vals)
            vals.to.match <- src.cdf(src)
        } 
        
        else if (match.method == "dens") {                      # Using kernel density
            kn.vals <- ref.cdf$eval.points;         fn.vals <- ref.cdf$estimate
            vals.to.match <- predict(src.cdf, x = src)
        } 
        
        else {
            stop("Invalid option for estimating distribution! Please use 'hist' for histogram or 'dens' for kernel density.")
        }
        
        ## Perform mapping...
        if (missing(lims))
            lims <- c(min(src), max(src))
        
        matched <- match.func(func.knots = kn.vals, func.vals = fn.vals, new.knots = vals.to.match, knot.lims = lims)
        matched
    }
    
    
    ## Predictive modeling...
    RF.predict <- function(x.train, y.train, x.test, y.lims, n.tree = 300, random.seed = NULL, ...) {
        if (!require(randomForest)) { library(randomForest) }       # Load package
        if (missing(y.lims)) { y.lims <- c(min(y.train), max(y.train)) }
        set.seed(seed = random.seed)                                # For reproducibility
        
        ## Define model & perform prediction...
        model  <- randomForest(x = x.train, y = y.train, ntree = n.tree, replace = TRUE, ...)
        y.pred <- conf.lims(predict(model, x.test), lims = y.lims)
        y.pred
    }
    
    
    ######## MAIN ##############################################################
    
    ## Define datasets...
    X1 <- norm.data(target.set[["X"]]);      y1 <- target.set[["y"]]
    X2 <- norm.data(source.set[["X"]]);      y2 <- source.set[["y"]]
    n.feat <- ncol(X1);                      data.lims <- c(0, 1)
    
    
    ## Distribution matching for predictors...
    X2.map <- lapply(1:n.feat, function(j) { 
        dist.match(X1[, j], ref = X2[, j], match.method = method, samp.size = size, lims = data.lims)
        })
    X2.map <- as.data.frame(X2.map);    dimnames(X2.map) <- dimnames(X1)
    # rownames(X2.map) <- rownames(X1);   colnames(X2.map) <- colnames(X1)
    
    
    ## Perform prediction & map back to original space...
    y2.pred.map <- RF.predict(x.train = X2, y.train = y2, x.test = X2.map, y.lims = data.lims, n.tree = 200, random.seed = seed)
    y2.dist <- get.cum.dist(y2, sample.size = size, dist.method = method, grid.size = 1e3)
    
    y1.pred <- dist.match(y2.pred.map, ref = y1, src.cdf = y2.dist, match.method = method, samp.size = size, lims = data.lims)
    y1.pred <- conf.lims(y1.pred, lims = data.lims);      names(y1.pred) <- names(y2.pred.map)
    
    
    ## Return output objects...
    if (pred.opt) {
        return( list("mapped" = y1.pred, "unmapped" = y2.pred.map) )
    } else {
        return( y1.pred )
    }

}