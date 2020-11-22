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


#### Functions...
printf <- function(..., end = "\n"){
  if ((nargs() > 1) & (grepl(list(...)[1], pattern = "%")))
    cat(sprintf(...), end)
  else
    cat(..., end)
}

norm01 <- function(x) (x - min(x)) / diff(range(x))
norm.data <- function(df) as.data.frame(apply(df, MARGIN = 2, norm01))

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
      stop("Invalid performance measure!")
  }
  
  perf.vals
}

# confine.in.lims <- function(y, lims = c(0, 1)) {
#   y[y < lims[1]] <- lims[1];      y[y > lims[2]] <- lims[2]
# }


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


#### Get data for a biomarker...
biomarkers <- colnames(Ydata1);       q <- length(biomarkers)


## Function for selecting top 'm' genes...
get.top.genes <- function(ranks, m_top = 150, print.opt = FALSE) {
  ## Get data for common top 'm' genes...
  nGN <- 300;     gnRank <- intersect(ranks[1:nGN, 2], ranks[1:nGN, 3])
  nI <- 0;        m <- length(gnRank);     m0 <- m
  while(m < m_top) {
    nI <- nI + 1;       nGN <- nGN + 100
    gnRank <- intersect(ranks[1:nGN, 2], ranks[1:nGN, 3]);    m <- length(gnRank)
  }
  
  ## Print results...
  if (print.opt)
    printf("#top genes chosen = %d (nGN = %d, nI = %d, m0 = %d)", m, nGN, nI, m0)
  
  ## Sort the ranks...
  gnRank <- sort(gnRank, decreasing = FALSE)
}
  

## Get results for all biomarkers...
source("dist.match.trans.learn.R")      ## Load function

run <- function(q.run, random.seed) {
# q.run <- 7                     # drug idx
# random.seed <- 4321            # 0, 654321, 4321

perf.mes <- c("NRMSE", "NMAE", "SCC")  
results.all <- list(data.frame("DMTL" = double(), "DMTL_SS" = double(), "BL" = double()), 
                    data.frame("DMTL" = double(), "DMTL_SS" = double(), "BL" = double()), 
                    data.frame("DMTL" = double(), "DMTL_SS" = double(), "BL" = double()), 
                    "genes" = data.frame("num.genes" = double()))
names(results.all)[1:3] <- perf.mes

for (k in q.run) {
  
  ## Select biomarker... 
  bmChosen <- biomarkers[k];      #printf("\nChosen biomarker = %s", bmChosen)
  ranks    <- cbind(rank1[, bmChosen], rank2[, bmChosen], rank3[, bmChosen])
  gnRank   <- get.top.genes(ranks, m_top = 150, print.opt = FALSE);      m <- length(gnRank)
  
  
  ## Prepare datasets...
  X1 <- Xdata1[, gnRank];               X2 <- rbind(Xdata2[, gnRank], Xdata3[, gnRank])
  Y1 <- norm01(Ydata1[, bmChosen]);     Y2 <- norm01(c(Ydata2[, bmChosen], Ydata3[, bmChosen]))
  
  
  ## DMTL model...
  prediction <- dist.match.trans.learn(target.set = list("X" = X1, "y" = Y1), source.set = list("X" = X2, "y" = Y2), 
                                       seed = random.seed, pred.opt = TRUE)
  Y1.pred <- prediction$mapped;     Y1.pred.src <- prediction$unmapped
  
  
  ## Baseline model...
  set.seed(random.seed)
  RF.base <- randomForest(x = norm.data(X2), y = Y2, ntree = 200, mtry = 5, replace = TRUE)
  Y1.pred.base <- predict(RF.base, norm.data(X1))
  Y1.pred.base[Y1.pred.base < 0] <- 0;      Y1.pred.base[Y1.pred.base > 1] <- 1
  
  
  ## Results df...
  # printf("After prediction: NRMSE = %0.4f, NMAE = %0.4f", NRMSE, NMAE)
  
  results <- data.frame("DMTL"    = calc.perf(Y1, Y1.pred, measures = perf.mes), 
                        "DMTL_SS" = calc.perf(Y1, Y1.pred.src, measures = perf.mes), 
                        "BL"      = calc.perf(Y1, Y1.pred.base, measures = perf.mes), row.names = perf.mes)
  # printf("Results = \n");   print(results)
  
  
  ## Save results in df...
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
printf("\nResults summary = ");    print(results.all$table)

results.all
}

# source("dist.match.trans.learn.R")      ## Load function
results.all <- run(q.run = 1:q, random.seed = 531)
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
