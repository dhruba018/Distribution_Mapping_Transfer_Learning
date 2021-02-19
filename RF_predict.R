## Perform predictive modeling using Random Forest regression 
## 
## Dependency: randomForest 
## Dependency_own: lambda_functions 
################################################################################

RF_predict <- function(x_train, y_train, x_test, y_lims, optimize = FALSE, n_tree = 300, m_try = 0.3333, random_seed = NULL, ...) {
    
    source("C:/Users/SRDhruba/Dropbox (Personal)/ResearchWork/Rtest/lambda_functions.R")
    
    ## Initital checks...
    # if (!require(randomForest))                 # Load package
    #     library(randomForest)
    
    if (missing(y_lims))
        y_lims <- c(min(y_train), max(y_train))
    
    if (m_try > 1 | m_try < 0)
        stop("Invalid value! Please choose a value between 0 and 1 (fraction of the features)!")
    
    
    ## Define model & perform prediction...
    set.seed(random_seed)                       # For reproducibility
    
    model_train <- function(x, y, method = "rf") {
        ctrl <- caret::trainControl(method = "cv", number = 10, search = "random")
        na <- which(is.na(y))
        if(length(na) > 0) {
            x <- as.matrix(x[-na, ]);    y <- y[-na]
        } else {
            x <- as.matrix(x)
        }
        # print(dim(x))
        Tune <- caret::train(x, y, method = method, tuneLength = 30, trControl = ctrl, preProc = NULL)
        # c("center", "scale", "medianImpute"))
        return(Tune)
    }
    
    if (optimize) {
        Forest <- model_train(x = x_train, y = y_train, method = "rf")
    } else {
        m_try  <- round(m_try * ncol(x_train))
        Forest <- randomForest::randomForest(x = x_train, y = y_train, ntree = n_tree, mtry = m_try, replace = TRUE, ...)
    }
    
    y_pred <- confined(stats::predict(Forest, x_test), lims = y_lims)
    y_pred
    
}