setwd("C:/Users/SRDhruba/Dropbox (Personal)/ResearchWork/Rtest/")

library(DMTL)
# library(CORElearn)
# library(randomForest)


### Read data...
load("./Data/ccle_gdsc_common_processed_SRD_19_Jan_2021.RData")
ccle_exp <- as.data.frame(ccle_exp, stringsAsFactors = FALSE)
ccle_auc <- as.data.frame(ccle_auc, stringsAsFactors = FALSE)
gdsc_exp <- as.data.frame(gdsc_exp, stringsAsFactors = FALSE)
gdsc_auc <- as.data.frame(gdsc_auc, stringsAsFactors = FALSE)


printf <- function(...) cat(sprintf(...), "\n")

## Files...
# source("dist_match_trans_learn.R");     source("dist_match.R");           source("get_dist_est.R")
# source("RF_predict.R");                 source("lambda_functions.R")


### Drug-wise analysis...
metrics <- c("NRMSE", "NMAE", "PCC", "SCC")
# results <- lapply(metrics, function(x) {
#   matrix(NA, nrow = length(drugs), ncol = 3, dimnames = list(drugs, c("DMTL", "DMTL_SS", "BL")))
# });    names(results) <- metrics


## Pick drug...
# drug <- drugs[1];
run <- function(drug, seed = NULL, K = 100, data_sw = 1, use_density = FALSE, verbose = TRUE) {

if (verbose)  printf("\nDrug = %s", drug)

ccle_idx <- which(!is.na(ccle_auc[, drug]));        ccle_cell <- ccle_cell_lines[ccle_idx]
ccle_set <- cbind(ccle_exp[ccle_idx, ], auc = ccle_auc[ccle_idx, drug])
gdsc_idx <- which(!is.na(gdsc_auc[, drug]));        gdsc_cell <- gdsc_cell_lines$cell[gdsc_idx]
gdsc_set <- cbind(gdsc_exp[gdsc_idx, ], auc = gdsc_auc[gdsc_idx, drug])
# print(ccle_set[1:8, c(1:8, ncol(ccle_set))]);       print(gdsc_set[1:8, c(1:8, ncol(gdsc_set))])


### Modeling...
## Switch for choosing datasets...
##  1 : {target = gdsc, source = ccle},   2 : {target = ccle, source = gdsc}
# data_sw <- 1

## Get sets.
load("./Data/relieff_fs.RData")
# K <- 100                                                        # No. of features

if (data_sw == 1) {                                             # target = gdsc, source = ccle
  if (verbose)  printf("Target = GDSC, Source = CCLE: Using top %d CCLE features", K)
  ccle_fs <- relieff_fs$ccle$ranks[, drug];       features <- ccle_fs[1:K]

  X1 <- gdsc_set[, features];                     Y1 <- gdsc_set$auc
  X2 <- ccle_set[, features];                     Y2 <- ccle_set$auc
} else {                                                        # target = ccle, source = gdsc
  if (verbose)  printf("Target = CCLE, Source = GDSC: Using top %d GDSC features", K)
  gdsc_fs <- relieff_fs$gdsc$ranks[, drug];       features <- gdsc_fs[1:K]

  X1 <- ccle_set[, features];                     Y1 <- ccle_set$auc
  X2 <- gdsc_set[, features];                     Y2 <- gdsc_set$auc
}


# seed <- 4

## DMTL.
prediction <- DMTL(target_set = list("X" = X1, "y" = Y1), source_set = list("X" = X2, "y" = Y2),
                   use_density = use_density, sample_size = 1e3, random_seed = seed, all_pred = TRUE)
Y1_pred <- prediction$target;     Y1_pred_src <- prediction$source
# prediction <- DMTL(target_set = list("X" = X1, "y" = Y1), source_set = list("X" = X2, "y" = Y2),
#                    method = "hist", size = 1e3, seed = seed, pred_all = TRUE)
# Y1_pred <- prediction$mapped;     Y1_pred_src <- prediction$unmapped


## Baseline.
Y1_pred_base <- RF_predict(x_train = norm_data(X2), y_train = Y2, x_test = norm_data(X1), lims = c(0, 1),
                           n_tree = 200, m_try = 0.4, seed = seed)
# Y1_pred_base <- RF_predict(x_train = norm_data(X2), y_train = Y2, x_test = norm_data(X1), y_lims = c(0, 1),
#                            n_tree = 200, m_try = 0.4, random_seed = seed)


## Performance evaluation.
perf_mat <- rbind(DMTL = performance(Y1, Y1_pred,      measures = metrics),
                  DMTL_SS = performance(Y1, Y1_pred_src,  measures = metrics),
                  BL   = performance(Y1, Y1_pred_base, measures = metrics))

if (verbose)  print(perf_mat)

perf_mat
}



### Run code...
run_all <- function(seed = NULL, K = 100, use_density = FALSE, data_sw = 1, progress = FALSE, verbose = FALSE) {
results <- lapply(metrics, function(x) {
  matrix(NA, nrow = length(drugs), ncol = 3, dimnames = list(drugs, c("DMTL", "DMTL_SS", "BL")))
});    names(results) <- metrics


if (progress)   pb <- progress::progress_bar$new(format = "[:bar] :percent eta: :eta", total = length(drugs), width = 64)
for (drug in drugs) {
if (progress)   pb$tick()
perf_mat <- run(drug = drug, seed = seed, K = K, use_density = use_density, data_sw = data_sw, verbose = verbose)

for (met in metrics) {
  results[[met]][drug, ] <- perf_mat[, met]
}
}

results <- lapply(metrics, function(x) rbind(results[[x]], mean = colMeans(results[[x]][drugs, ])));   names(results) <- metrics
res_sum <- Reduce(lapply(metrics, function(x) results[[x]]["mean", ]), f = rbind);    rownames(res_sum) <- metrics
print(res_sum)

# run(drug = drugs[1], seed = 97531)

print( c(sum(results$NRMSE[drugs, "DMTL"] > 1), sum(results$NMAE[drugs, "DMTL"] > 1),
         sum(results$PCC[drugs, "DMTL"] < 0.2), sum(results$SCC[drugs, "DMTL"] < 0.2)) )

results
}

# results <- run_all(seed = 4, K = 150)
results <- run_all(seed = 4, K = 200, use_density = FALSE, progress = TRUE)
results <- run_all(seed = 4, K = 200, use_density = TRUE, progress = TRUE)

results2 <- run_all(seed = 4, K = 200, use_density = FALSE, data_sw = 2, progress = TRUE, verbose = TRUE)
results2 <- run_all(seed = 4, K = 200, use_density = TRUE, data_sw = 2, progress = TRUE, verbose = TRUE)


