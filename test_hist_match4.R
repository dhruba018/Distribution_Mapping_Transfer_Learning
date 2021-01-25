# rm(list = ls())

## Set-up system path... 
info <- Sys.getenv(c("USERNAME", "HOMEPATH"))
if (info["USERNAME"] == "SRDhruba"){
  info["DIRPATH"] <- sprintf("%s\\Dropbox%s\\ResearchWork\\Rtest\\", info["HOMEPATH"], " (Personal)")
} else {
  info["DIRPATH"] <- sprintf("%s\\Dropbox%s\\ResearchWork\\Rtest\\", info["HOMEPATH"], "")
}
setwd(info["DIRPATH"]);       cat("Current system path = ", getwd(), "\n")


## Packages...
library(ggplot2)
library(ggpubr)
library(randomForest)
library(ks)
library(progress)


#### Functions...
printf <- function(..., end = "\n") {
  if ((nargs() > 1) & (grepl(list(...)[1], pattern = "%")))
    cat(sprintf(...), end)
  else
    cat(..., end)
}

norm01    <- function(z) { z <- if (min(z)) z - min(z) else z;   z <- z / max(z);  z }
norm.data <- function(df) as.data.frame(apply(df, MARGIN = 2, norm01))


### Pick top genes...
get.top.genes <- function(ranks, m.top = 150, verbose = FALSE) {
  
  ## Initialization...
  nI <- 0;        nGN <- 300
  gene.rank <- intersect(ranks[1:nGN, 1], ranks[1:nGN, 2])
  m  <- length(gene.rank);     m0 <- if (verbose) m
  
  ## Run iterations...
  while(m < m.top) {
    nI <- nI + 1;       nGN <- nGN + 100
    gene.rank <- intersect(ranks[1:nGN, 1], ranks[1:nGN, 2])
    m  <- length(gene.rank)
  }
  gene.rank <- sort(gene.rank, decreasing = FALSE)    # Sort ranks
  
  ## Print results...
  if (verbose)
    printf("#top genes chosen = %d (nGN = %d, nI = %d, m0 = %d)", m, nGN, nI, m0)
  
  gene.rank
}


### Get model performance...
calc.perf <- function(y, y.pred, measures = c("NRMSE", "NMAE", "SCC")) {
  
  ## Initialize results array...
  perf.vals <- c("sq.err" = NA, "abs.err" = NA, "cor.coef" = NA)
  
  for (mes in measures) {
    
    ## Calculate squared error...
    if (grepl(pattern = "MSE", mes, ignore.case = TRUE)) {
      num <- mean((y - y.pred)^2)
      den <- if (mes == "NRMSE") mean((y - mean(y))^2) else 1
      pow <- if (mes == "MSE") 1 else 0.5
      perf.vals["sq.err"] <- (num / den)^pow
    } 
    
    ## Calculate absolute error...
    else if (grepl(pattern = "MAE", mes, ignore.case = TRUE)) {
      num <- mean(abs(y - y.pred))
      den <- if (mes == "NMAE") mean(abs(y - mean(y))) else 1
      perf.vals["abs.err"] <- num / den
    } 
    
    ## Calculate similarity measures...
    else if (grepl(pattern = "CC", mes, ignore.case = TRUE)) {
      alg <- if (mes == "SCC") "spearman" else "pearson"
      perf.vals["cor.coef"] <- cor(y, y.pred, method = alg)
    } 
    
    ## Doesn't match any...
    else 
      stop("Invalid performance measure! Please use common variants of MSE, MAE or CC (correlation coefficient).")
  }
  
  perf.vals
}


#### Read tumor-cell line data...
Xdata1 <- read.table("Data/BRCA_gene_expression_METABRIC_26_Oct_2020.txt", sep = "\t", header = TRUE)
Xdata2 <- read.table("Data/BRCA_gene_expression_CCLE_26_Oct_2020.txt", sep = "\t", header = TRUE)
Xdata3 <- read.table("Data/BRCA_gene_expression_GDSC_26_Oct_2020.txt", sep = "\t", header = TRUE)

Ydata1 <- read.table("Data/BRCA_biomarker_expression_METABRIC_26_Oct_2020.txt", sep = "\t", header = TRUE)
Ydata2 <- read.table("Data/BRCA_biomarker_expression_CCLE_26_Oct_2020.txt", sep = "\t", header = TRUE)
Ydata3 <- read.table("Data/BRCA_biomarker_expression_GDSC_26_Oct_2020.txt", sep = "\t", header = TRUE)

rank1 <- read.table("Data/BRCA_biomarker_ranks_METABRIC_27_Oct_2020.txt", sep = "\t", header = TRUE)
rank2 <- read.table("Data/BRCA_biomarker_ranks_CCLE_27_Oct_2020.txt", sep = "\t", header = TRUE)
rank3 <- read.table("Data/BRCA_biomarker_ranks_GDSC_27_Oct_2020.txt", sep = "\t", header = TRUE)

biomarkers <- colnames(Ydata1);       q <- length(biomarkers)


## Get results for all biomarkers...
# source("dist.match.trans.learn.R")      ## Load function
# source("dist_match_trans_learn.R")      ## Load function

run <- function(q.run, n.feat, random.seed, method.opt) {
# q.run <- 1                     # drug idx
# random.seed <- 4321            # 0, 654321, 4321
# method.opt <- "dens"           # hist, dens

source("RF_predict.R")          # Random forest modeling

perf.mes <- c("NRMSE", "NMAE", "SCC")  
results.all <- list(data.frame("DMTL" = double(), "DMTL_SS" = double(), "BL" = double()), 
                    data.frame("DMTL" = double(), "DMTL_SS" = double(), "BL" = double()), 
                    data.frame("DMTL" = double(), "DMTL_SS" = double(), "BL" = double()), 
                    "genes" = data.frame("num.genes" = double()))
names(results.all)[1:3] <- perf.mes

pb <- progress_bar$new(format = "  running [:bar] :percent eta: :eta", total = length(q.run), clear = FALSE, width = 64)
pb$tick(0)

for (k in q.run) {
  
  pb$tick()
  
  ## Select biomarker... 
  bmChosen <- biomarkers[k];      #printf("\nChosen biomarker = %s", bmChosen)
  ranks    <- cbind(rank1[, bmChosen], rank2[, bmChosen], rank3[, bmChosen])
  gnRank   <- get.top.genes(ranks[, 2:3], m.top = n.feat, verbose = FALSE);      m <- length(gnRank)
  
  
  ## Prepare datasets...
  X1 <- Xdata1[, gnRank];               X2 <- rbind(Xdata2[, gnRank], Xdata3[, gnRank])
  Y1 <- norm01(Ydata1[, bmChosen]);     Y2 <- norm01(c(Ydata2[, bmChosen], Ydata3[, bmChosen]))
  
  
  ## DMTL model...
  prediction <- DMTL(target_set = list("X" = X1, "y" = Y1), source_set = list("X" = X2, "y" = Y2), 
                     method = method.opt, seed = random.seed, pred_all = TRUE)
  Y1.pred <- prediction$mapped;     Y1.pred.src <- prediction$unmapped
  
  
  ## Baseline model...
  # set.seed(random.seed)
  # RF.base <- randomForest(x = norm.data(X2), y = Y2, ntree = 200, mtry = 5, replace = TRUE)
  # Y1.pred.base <- predict(RF.base, norm.data(X1))
  # Y1.pred.base[Y1.pred.base < 0] <- 0;      Y1.pred.base[Y1.pred.base > 1] <- 1
  # 
  Y1.pred.base <- RF_predict(x_train = norm.data(X2), y_train = Y2, x_test = norm.data(X1), 
                             n_tree = 200, m_try = 0.4, random_seed = random.seed)
  
  
  ## Generate & save results...
  results <- data.frame("DMTL"    = calc.perf(Y1, Y1.pred, measures = perf.mes), 
                        "DMTL_SS" = calc.perf(Y1, Y1.pred.src, measures = perf.mes), 
                        "BL"      = calc.perf(Y1, Y1.pred.base, measures = perf.mes), row.names = perf.mes)
  
  ## Print option...
  if (length(q.run) == 1) { printf("\nResults for %s using top %d features = ", bmChosen, n.feat);     print(results) }
  
  results.all[[perf.mes[1]]][bmChosen, ] <- results[perf.mes[1], ]
  results.all[[perf.mes[2]]][bmChosen, ] <- results[perf.mes[2], ]
  results.all[[perf.mes[3]]][bmChosen, ] <- results[perf.mes[3], ]
  results.all$genes[bmChosen, ] <- m
}

## Calculate mean performance...
results.all[[perf.mes[1]]]["Mean", ] <- colMeans(results.all[[perf.mes[1]]][biomarkers, ], na.rm = TRUE)
results.all[[perf.mes[2]]]["Mean", ] <- colMeans(results.all[[perf.mes[2]]][biomarkers, ], na.rm = TRUE)
results.all[[perf.mes[3]]]["Mean", ] <- colMeans(results.all[[perf.mes[3]]][biomarkers, ], na.rm = TRUE)
results.all$genes["Mean", ]          <- mean(results.all$genes[biomarkers, ], na.rm = TRUE)

results.all[["table"]] <- rbind(results.all[[perf.mes[1]]]["Mean", ], results.all[[perf.mes[2]]]["Mean", ], 
                                results.all[[perf.mes[3]]]["Mean", ])
rownames(results.all$table) <- perf.mes

## Print options...
if (length(q.run) > 1) { printf("\nResults summary for top %d features = ", n.feat);    print(results.all$table) }

results.all
}

# source("dist.match.trans.learn.R")      ## Load function
source("dist_match_trans_learn.R")      ## Load function
results.all <- run(q.run = 1:q, n.feat = 50, random.seed = 97531, method.opt = "hist")
c(sum(results.all$NRMSE$DMTL >= 1), sum(results.all$NMAE$DMTL >= 1), sum(abs(results.all$SCC$DMTL) <= 0.2))


# ## Write in temporary file...
# write.in.file <- function() {
# write.table(results.all$NRMSE, file = sprintf("results_temp_%s.csv", format(Sys.Date(), "%d_%b_%Y")),
#             sep = "\t", row.names = TRUE, col.names = TRUE)
# write.table(results.all$NMAE, file = sprintf("results_temp_%s.csv", format(Sys.Date(), "%d_%b_%Y")),
#             sep = "\t", append = TRUE, row.names = TRUE, col.names = TRUE)
# write.table(results.all$SCC, file = sprintf("results_temp_%s.csv", format(Sys.Date(), "%d_%b_%Y")),
#             sep = "\t", append = TRUE, row.names = TRUE, col.names = TRUE)
# }
# write.in.file()



# ##
MB.df <- as.data.frame(cbind(Y1, Y1.pred, Y1.pred.src, Y1.pred.base), row.names = rownames(Ydata1))
CG.df <- as.data.frame(Y2, row.names = c(rownames(Ydata2), rownames(Ydata3)))

plot.dist <- function() {
gg.pp <- list()
gg.pp[["Tumor"]] <- ggplot(MB.df, aes(x = Y1)) + geom_density(fill = "cyan4", alpha = 0.7) + 
                      geom_histogram(aes(y = ..density..), bins = 100, color = "gray2", 
                                     fill = "brown2", alpha = 0.3) + theme_light() + xlab("") + 
                      ylab("") + ggtitle("Tumor") + theme(plot.title = element_text(hjust = 0.5))
gg.pp[["Cell line"]] <- ggplot(CG.df, aes(x = Y2)) + geom_density(fill = "brown4", alpha = 0.7) + 
                          geom_histogram(aes(y = ..density..), bins = 100, color = "gray2", 
                                         fill = "cyan2", alpha = 0.3) + theme_light() + xlab("") + 
                          ylab("") + ggtitle("Cell line") + theme(plot.title = element_text(hjust = 0.5))
gg.pp[["Tumor.Pred"]] <- ggplot(MB.df, aes(x = Y1.pred)) + geom_density(fill = "cyan4", alpha = 0.7) + 
                          geom_histogram(aes(y = ..density..), bins = 100, color = "gray2", 
                                         fill = "brown2", alpha = 0.3) + theme_light() + xlab("") + 
                          ylab("") + ggtitle("Tumor (Predicted)") + theme(plot.title = element_text(hjust = 0.5))
gg.pp[["Tumor.Pred.Src"]] <- ggplot(MB.df, aes(x = Y1.pred.src)) + geom_density(fill = "brown4", alpha = 0.7) + 
                              geom_histogram(aes(y = ..density..), bins = 100, color = "gray2", 
                                             fill = "cyan2", alpha = 0.3) + theme_light() + xlab("") + 
                              ylab("") + ggtitle("Tumor (Predicted - Source)") + 
                              theme(plot.title = element_text(hjust = 0.5))

gg.pp[c("ncol", "nrow")] <- list(2, 2)
annotate_figure(do.call(ggarrange, gg.pp), top = text_grob(bmChosen, face = "bold"))
}

plot.dist()
