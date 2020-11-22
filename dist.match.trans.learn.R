dist.match.trans.learn <- function(target.set, source.set, method = "hist", size = 1e3, seed = NULL, pred.opt = FALSE){
    
    ## Initial check...
    if (ncol(target.set[["X"]]) != ncol(source.set[["X"]]))
        stop("Source and target set covariates must have the same number of features!")
    
    if (nrow(source.set[["X"]]) != length(source.set[["y"]]))
        stop("Source sets must have the same number of samples!")
    
    out01 <- function(y) (min(y) < 0) & (max(y) > 1)
    if (out01(target.set[["y"]]) | out01(source.set[["y"]]))
        stop("Response data must be normalized in [0, 1].")
    
    
    ## Functions...
    norm01    <- function(z) {z <- if (min(z)) z - min(z);   z <- z / max(z);     z }
    norm.data <- function(df) as.data.frame(apply(df, MARGIN = 2, norm01))
    zscore    <- function(df) as.data.frame(apply(df, MARGIN = 2, scale))
    conf.lims <- function(y, lims = c(0, 1)){ y[y < lims[1]] <- lims[1];    y[y > lims[2]] <- lims[2];   y}
    
    ## Estimate distributions...
    get.cum.dist <- function(x, sample.size = 1e6, dist.method = "hist") {
        dist.method = tolower(dist.method)
        if (dist.method == "hist") {
            xx <- sample(x, size = sample.size, replace = TRUE)
            x.cdf <- ecdf(xx)
        } else if (dist.method == "kde") {
            # x.cdf <- kCDF
        }
    }
    
    ## Distribution matching...
    dist.match <- function(src, ref, src.cdf, ref.cdf, lims, match.method = "hist", samp.size = 1e6) {
        ## Get distribution...
        if (missing(src.cdf))
            src.cdf <- get.cum.dist(src, sample.size = samp.size, dist.method = match.method)

        if (missing(ref.cdf))
            ref.cdf <- get.cum.dist(ref, sample.size = samp.size, dist.method = match.method)

        ## Get mapping...
        if (missing(lims))
            lims <- c(floor(min(src)), ceiling(max(src)))
        
        kn <- knots(ref.cdf);       fn <- ref.cdf(kn)
        map <- approxfun(x = fn, y = kn, method = "linear", yleft = lims[1], yright = lims[2], ties = "ordered")
        
        ## Mapped data...
        matched <- map(src.cdf(src))
    }
    
    ## Predictive modeling...
    RF.prediction <- function(x.train, y.train, x.test, n.tree = 200, random.seed = NULL, ...) {
        ## Set up model...
        if (!require(randomForest))             # Load package
            library(randomForest)
        set.seed(seed = random.seed)            # For reproducibility
        
        ## Define model & perform prediction...
        forest <- randomForest(x = x.train, y = y.train, ntree = n.tree, replace = TRUE, ...)
        y.pred <- conf.lims(predict(forest, x.test), lims = c(0, 1))
        y.pred
    }
    
    
    ######## MAIN ##############################################################
    
    ## Define datasets...
    X1 <- norm.data(target.set[["X"]]);      y1 <- target.set[["y"]]
    X2 <- norm.data(source.set[["X"]]);      y2 <- source.set[["y"]]
    n.feat <- ncol(X1)
    
    ## Distribution matching for predictors...
    X2.map <- lapply(1:n.feat, function(j) { 
        dist.match(X1[, j], ref = X2[, j], match.method = method, samp.size = size)
        })
    X2.map <- as.data.frame(X2.map, col.names = colnames(X1), row.names = rownames(X1))
    
    ## Perform prediction & map back to original space...
    y2.pred.map <- RF.prediction(x.train = X2, y.train = y2, x.test = X2.map, n.tree = 200, random.seed = seed)
    
    y1.pred <- dist.match(y2.pred.map, ref = y1, src.cdf = get.cum.dist(y2, sample.size = size, dist.method = method), 
                          match.method = method, samp.size = size)
    y1.pred <- conf.lims(y1.pred, lims = c(0, 1));      names(y1.pred) <- names(y2.pred.map)
    
    ## Return output objects...
    if (pred.opt) {
        return( list("mapped" = y1.pred, "unmapped" = y2.pred.map) )
    } else {
        return( y1.pred )
    }

}