# rm(list = ls())

## Set-up system path...
info <- Sys.getenv(c("USERNAME", "HOMEPATH"))
if (info["USERNAME"] == "SRDhruba"){
  info["DIRPATH"] <- sprintf("%s\\Dropbox (Personal)\\ResearchWork\\Rtest\\", info["HOMEPATH"])
} else {
  info["DIRPATH"] <- sprintf("%s\\Dropbox\\ResearchWork\\Rtest\\", info["HOMEPATH"])
}
setwd(info["DIRPATH"]);       cat("Current system path = ", getwd(), "\n")


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
  perf.vals <- c()

  for (mes in measures) {

    ## Calculate squared error...
    if (grepl(pattern = "MSE", mes, ignore.case = TRUE)) {
      num <- mean((y - y.pred)^2)
      den <- if (mes == "NRMSE") mean((y - mean(y))^2) else 1
      pow <- if (mes == "MSE") 1 else 0.5
      perf.vals[mes] <- (num / den)^pow
    }

    ## Calculate absolute error...
    else if (grepl(pattern = "MAE", mes, ignore.case = TRUE)) {
      num <- mean(abs(y - y.pred))
      den <- if (mes == "NMAE") mean(abs(y - mean(y))) else 1
      perf.vals[mes] <- num / den
    }

    ## Calculate similarity measures...
    else if (grepl(pattern = "CC", mes, ignore.case = TRUE)) {
      alg <- if (mes == "SCC") "spearman" else "pearson"
      perf.vals[mes] <- cor(y, y.pred, method = alg)
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
run <- function(q.run, n.feat, random.seed, density.opt = FALSE, model = "RF", optimize = FALSE) {
  perf.mes <- c("NRMSE", "NMAE", "PCC", "SCC")
  results.all <- lapply(perf.mes, function(mes) data.frame("DMTL" = double(), "DMTL_SS" = double(), "BL" = double()))
  names(results.all) <- perf.mes;       results.all[["genes"]] <- data.frame("num.genes" = double())


  ## Prints progress...
  pb <- progress_bar$new(format = "  running [:bar] :percent eta: :eta", total = length(q.run), clear = FALSE, width = 64)
  pb$tick(0)

  for (k in q.run) {

    pb$tick()

    # ## When testing for a single drug...
    # k <- 1;     n.feat <- 150;    density.opt <- FALSE;     random.seed <- 7531

    ## Select biomarker...
    bmChosen <- biomarkers[k];      #printf("\nChosen biomarker = %s", bmChosen)
    ranks    <- cbind(rank1[, bmChosen], rank2[, bmChosen], rank3[, bmChosen])
    gnRank   <- get.top.genes(ranks[, 2:3], m.top = n.feat, verbose = FALSE);      m <- length(gnRank)


    ## Prepare datasets...
    X1 <- Xdata1[, gnRank];               X2 <- rbind(Xdata2[, gnRank], Xdata3[, gnRank])
    Y1 <- norm01(Ydata1[, bmChosen]);     Y2 <- norm01(c(Ydata2[, bmChosen], Ydata3[, bmChosen]))
    names(Y1) <- rownames(X1);            names(Y2) <- rownames(X2)


    ## DMTL model...
    prediction <- DMTL(target_set = list("X" = X1, "y" = Y1), source_set = list("X" = X2, "y" = Y2), pred_model = model,
                       model_optimize = optimize, use_density = density.opt, random_seed = random.seed, all_pred = TRUE)
    Y1.pred <- prediction$target;     Y1.pred.src <- prediction$source


    ## Baseline model...
    Y1.pred.base <- if (model == "RF") {
      RF_predict(x_train = norm.data(X2), y_train = Y2, x_test = norm.data(X1), lims = c(0, 1), optimize = optimize,
                 n_tree = 200, m_try = 0.4, seed = random.seed)
    } else if (model == "SVM") {
      SVM_predict(x_train = norm.data(X2), y_train = Y2, x_test = norm.data(X1), lims = c(0, 1), optimize = optimize,
                  kernel = "rbf", C = 2, eps = 0.01, kpar = list(sigma = 0.1), seed = random.seed)
    } else if (model == "EN") {
      EN_predict(x_train = norm.data(X2), y_train = Y2, x_test = norm.data(X1), lims = c(0, 1), optimize = optimize,
                 alpha = 0.8, seed = random.seed)
    }


    ## Generate & save results...
    results <- data.frame("DMTL"    = calc.perf(Y1, Y1.pred, measures = perf.mes),
                          "DMTL_SS" = calc.perf(Y1, Y1.pred.src, measures = perf.mes),
                          "BL"      = calc.perf(Y1, Y1.pred.base, measures = perf.mes), row.names = perf.mes)

    ## Print option...
    if (length(q.run) == 1) { printf("\nResults for %s using top %d features = ", bmChosen, n.feat);     print(results) }

    for (mes in perf.mes) { results.all[[mes]][bmChosen, ] <- results[mes, ] }
    results.all$genes[bmChosen, ] <- m
  }

  ## Calculate mean performance...
  for (mes in perf.mes) { results.all[[mes]]["Mean", ] <- colMeans(results.all[[mes]][biomarkers, ], na.rm = TRUE) }
  results.all$genes["Mean", ] <- mean(results.all$genes[biomarkers, ], na.rm = TRUE)

  results.all[["table"]] <- do.call(rbind, lapply(perf.mes, function(mes) results.all[[mes]]["Mean", ]))
  rownames(results.all$table) <- perf.mes

  ## Print options...
  if (length(q.run) > 1) { printf("\nResults summary for top %d features = ", n.feat);    print(results.all$table) }

  results.all
}

## Try out different models...
results.all.rf  <- run(q.run = 1:q, n.feat = 150, random.seed = 7531, density.opt = FALSE, model = "RF",  optimize = FALSE)
results.all.svm <- run(q.run = 1:q, n.feat = 150, random.seed = 7531, density.opt = FALSE, model = "SVM", optimize = FALSE)
results.all.en  <- run(q.run = 1:q, n.feat = 150, random.seed = 7531, density.opt = FALSE, model = "EN",  optimize = FALSE)


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


