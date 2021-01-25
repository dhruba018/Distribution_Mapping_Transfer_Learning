# rm(list = ls())

## Set-up system path...
PATH <- if (Sys.getenv("USERNAME") == "SRDhruba") {
  "\\Users\\SRDhruba\\Dropbox (Personal)\\ResearchWork\\Rtest\\"
} else {
  sprintf("%s\\Dropbox\\ResearchWork\\Rtest\\", Sys.getenv("HOMEPATH"))
}
setwd(PATH);       cat("Current system path = ", getwd(), "\n")


## Packages...
library(progress)
library(DMTL)


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
Xdata1 <- read.table("Data/LUAD_gene_expression_TCGA_06_Dec_2020.txt", sep = "\t", header = TRUE)
Xdata2 <- read.table("Data/LUSC_gene_expression_TCGA_06_Dec_2020.txt", sep = "\t", header = TRUE)
Xdata3 <- read.table("Data/NSCLC_gene_expression_CCLE_06_Dec_2020.txt", sep = "\t", header = TRUE)
Xdata4 <- read.table("Data/NSCLC_gene_expression_GDSC_06_Dec_2020.txt", sep = "\t", header = TRUE)

Ydata1 <- read.table("Data/LUAD_biomarker_expression_TCGA_06_Dec_2020.txt", sep = "\t", header = TRUE)
Ydata2 <- read.table("Data/LUSC_biomarker_expression_TCGA_06_Dec_2020.txt", sep = "\t", header = TRUE)
Ydata3 <- read.table("Data/NSCLC_biomarker_expression_CCLE_06_Dec_2020.txt", sep = "\t", header = TRUE)
Ydata4 <- read.table("Data/NSCLC_biomarker_expression_GDSC_06_Dec_2020.txt", sep = "\t", header = TRUE)

rank1 <- read.table("Data/LUAD_biomarker_ranks_TCGA_06_Dec_2020.txt", sep = "\t", header = TRUE)
rank2 <- read.table("Data/LUSC_biomarker_ranks_TCGA_06_Dec_2020.txt", sep = "\t", header = TRUE)
rank3 <- read.table("Data/NSCLC_biomarker_ranks_CCLE_06_Dec_2020.txt", sep = "\t", header = TRUE)
rank4 <- read.table("Data/NSCLC_biomarker_ranks_GDSC_06_Dec_2020.txt", sep = "\t", header = TRUE)

biomarkers <- colnames(Ydata1);       q <- length(biomarkers)


## Get results for all biomarkers...
# source("dist_match_trans_learn.R")      ## Load function

run <- function(q.run, n.feat, random.seed, density.opt) {
  # q.run <- 1                     # drug idx
  # random.seed <- 4321            # 0, 654321, 4321
  # method.opt <- "dens"           # hist, dens

  ## Save performance measures...
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
    ranks    <- cbind(rank1[, bmChosen], rank2[, bmChosen], rank3[, bmChosen], rank4[, bmChosen])
    gnRank   <- get.top.genes(ranks[, 3:4], m.top = n.feat, verbose = FALSE);      m <- length(gnRank)


    ## Prepare datasets...
    X1 <- rbind(Xdata1[, gnRank], Xdata2[, gnRank]);              X2 <- rbind(Xdata3[, gnRank], Xdata4[, gnRank])
    Y1 <- norm01(c(Ydata1[, bmChosen], Ydata2[, bmChosen]));      Y2 <- norm01(c(Ydata3[, bmChosen], Ydata4[, bmChosen]))


    ## DMTL model...
    prediction <- DMTL(target_set = list("X" = X1, "y" = Y1), source_set = list("X" = X2, "y" = Y2),
                       use_density = density.opt, random_seed = random.seed, all_pred = TRUE)
    Y1.pred <- prediction$target;     Y1.pred.src <- prediction$source


    ## Baseline model...
    # source("RF_predict.R")          # Random forest modeling

    Y1.pred.base <- RF_predict(x_train = norm.data(X2), y_train = Y2, x_test = norm.data(X1),
                               lims = c(0, 1), n_tree = 200, m_try = 0.4, random_seed = random.seed)


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

# source("dist_match_trans_learn.R")      ## Load function
results.all <- run(q.run = 1:q, n.feat = 50, random.seed = 97531, density.opt = FALSE)
c(sum(results.all$NRMSE$DMTL >= 1), sum(results.all$NMAE$DMTL >= 1), sum(abs(results.all$SCC$DMTL) <= 0.2))



