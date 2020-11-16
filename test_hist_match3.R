# rm(list = ls())

## Set-up system path... 
info <- Sys.getenv(c("USERNAME", "HOMEPATH"))
if (info["USERNAME"] == "SRDhruba"){
  info["DIRPATH"] <- sprintf("%s\\Dropbox%s\\ResearchWork\\Rtest\\", info["HOMEPATH"], " (Personal)")
} else {
  info["DIRPATH"] <- sprintf("%s\\Dropbox%s\\ResearchWork\\Rtest\\", info["HOMEPATH"], "")
}
setwd(info["DIRPATH"]);       cat("Current system path = ", getwd())


## Packages...
library(ggplot2)
library(ggpubr)
library(randomForest)


#### Functions...
printf <- function(..., end = "\n"){
  if ((nargs() > 1) & (grepl(list(... )[1], pattern = "%")))
    cat(sprintf(...), end)
  else
    cat(..., end)
}

dapply <- function(df, ...) as.data.frame(apply(df, ...))

norm01 <- function(x) (x - min(x)) / diff(range(x))

calc.err <- function(y, y.pred, measure = "MSE"){
  measure <- toupper(measure)
  if (grepl(pattern = "MSE", measure)){
    err <- mean((y - y.pred)^2)
    if (measure == "RMSE"){
      err <- sqrt(err)
    } else if (measure == "NRMSE"){
      err <- sqrt(err / mean((y - mean(y))^2))
    }
  } else if (grepl(pattern = "MAE", measure)){
    err <- mean(abs(y - y.pred))
    if (measure == "NMAE"){
      err <- err / mean(abs(y - mean(y)))
    }
  }
  err
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


#### Get data for a biomarker...
biomarkers <- colnames(Ydata1);       bmChosen <- biomarkers[6]
ranks <- cbind(rank1[, bmChosen], rank2[, bmChosen], rank3[, bmChosen])
printf("Chosen biomarker = %s", bmChosen)

## Get data for common top 'm' genes...
m_opt <- 150;       nGN <- 300;               gnRank <- intersect(ranks[1:nGN, 2], ranks[1:nGN, 3])
nI <- 0;            m0 <- length(gnRank);     m <- m0
while(m < m_opt){
  nI <- nI + 1;     nGN <- nGN + 100
  gnRank <- intersect(ranks[1:nGN, 2], ranks[1:nGN, 3]);      m <- length(gnRank)
}
gnRank <- sort(gnRank, decreasing = FALSE)
printf("#top genes chosen = %d (nGN = %d, nI = %d, m0 = %d)", m, nGN, nI, m0)


X1 <- Xdata1[, gnRank];               X2 <- rbind(Xdata2[, gnRank], Xdata3[, gnRank])
Y1 <- norm01(Ydata1[, bmChosen]);     Y2 <- norm01(c(Ydata2[, bmChosen], Ydata3[, bmChosen]))


## DMTL model...
source("dist.match.trans.learn.R")
Y1.pred <- dist.match.trans.learn(target.set = list("X" = X1, "y" = Y1), 
                                  source.set = list("X" = X2, "y" = Y2))


## Baseline model...
# set.seed(0)
RF.base <- randomForest(x = dapply(X2, MARGIN = 2, norm01), y = Y2, ntree = 200)
Y1.pred.base <- predict(RF.base, dapply(X1, MARGIN = 2, norm01))
Y1.pred.base[Y1.pred.base < 0] <- 0;    Y1.pred.base[Y1.pred.base > 1] <- 1


## Results df...
# printf("After prediction: NRMSE = %0.4f, NMAE = %0.4f", NRMSE, NMAE)

results <- data.frame("DMTL" = c(calc.err(Y1, Y1.pred, measure = "NRMSE"), calc.err(Y1, Y1.pred, measure = "NMAE"), 
                                 cor(Y1, Y1.pred, method = "spearman")), 
                      'BL' = c(calc.err(Y1, Y1.pred.base, measure = "NRMSE"), calc.err(Y1, Y1.pred.base, measure = "NMAE"), 
                               cor(Y1, Y1.pred.base, method = "spearman")), 
                      row.names = c("NRMSE", "NMAE", "SCC"))

printf("Results = \n");   print(results)


# ##
# gg.pp <- list()
# gg.pp[["Tumor"]] <- ggplot(MB.df, aes_string(x = bmChosen)) + geom_density(fill = "cyan4", alpha = 0.7) + theme_light() +
#                       xlab("") + ylab("") + ggtitle("Tumor") + theme(plot.title = element_text(hjust = 0.5))
# gg.pp[["Cell line"]] <- ggplot(CG.df, aes_string(x = bmChosen)) + geom_density(fill = "brown4", alpha = 0.7) + theme_light() +
#                           xlab("") + ylab("") + ggtitle("Cell line") + theme(plot.title = element_text(hjust = 0.5))
# gg.pp[["ncol"]] <- 2
# annotate_figure(do.call(ggarrange, gg.pp), top = text_grob(bmChosen, face = "bold"))
