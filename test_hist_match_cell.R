setwd("C:/Users/SRDhruba/Dropbox (Personal)/ResearchWork/Rtest/")

library(CORElearn)
# library(FSelector)
# library(randomForest)
library(DMTL)


### Read data...
load("./Data/ccle_gdsc_common_processed_SRD_19_Jan_2021.RData")
ccle_exp <- as.data.frame(ccle_exp, stringsAsFactors = FALSE)
ccle_auc <- as.data.frame(ccle_auc, stringsAsFactors = FALSE)
gdsc_exp <- as.data.frame(gdsc_exp, stringsAsFactors = FALSE)
gdsc_auc <- as.data.frame(gdsc_auc, stringsAsFactors = FALSE)


### Drug-wise analysis...
# relieff_fs <- 
#   list(ccle = list(weights = matrix(nrow = length(genes), ncol = length(drugs), dimnames = list(genes, drugs)), 
#                    ranks   = matrix(nrow = length(genes), ncol = length(drugs), dimnames = list(1:length(genes), drugs))), 
#        gdsc = list(weights = matrix(nrow = length(genes), ncol = length(drugs), dimnames = list(genes, drugs)), 
#                    ranks   = matrix(nrow = length(genes), ncol = length(drugs), dimnames = list(1:length(genes), drugs))))


## Files & functions...
# source("dist_match_trans_learn.R");     source("dist_match.R");           source("get_dist_est.R")
# source("RF_predict.R");                 source("lambda_functions.R")

COR   <- function(y_act, y_pred, sw = 1) cor(y_act, y_pred, method = if (sw == 1) "spearman" else method = "pearson")
NRMSE <- function(y_act, y_pred, sw = 1) { 
  sqrt(mean((y_act - y_pred)^2)) / (if (sw == 1) sqrt(mean((y_act - mean(y_act))^2)) else diff(range(y_act)))
}
NMAE  <- function(y_act, y_pred, sw = 1) {
  mean(abs(y_act - y_pred)) / (if (sw == 1) mean(abs(y_act - mean(y_act))) else diff(range(y_act)))
}



## Pick drug...
load("./Data/relieff_fs.RData")

# drug <- drugs[4];     

# run <- function(seed = NULL, K = 100, verbose = TRUE) {
# results_ccle <- matrix(nrow = length(drugs), ncol = 6, dimnames = 
#                          list(drugs, c("CC1", "NRMSE1", "NMAE1", "CC2", "NRMSE2", "NMAE2")))
# results_ccle <- list(DMTL = results_ccle, DMTL_SS = results_ccle, BL = results_ccle)
results_gdsc <- matrix(nrow = length(drugs), ncol = 6, dimnames = 
                         list(drugs, c("CC1", "NRMSE1", "NMAE1", "CC2", "NRMSE2", "NMAE2")))
results_gdsc <- list(DMTL = results_gdsc, DMTL_SS = results_gdsc, BL = results_gdsc)

# for (drug in drugs) {
# if (verbose) {
# print(sprintf("Drug = %s", drug))
# }

ccle_idx <- which(!is.na(ccle_auc[, drug]));        ccle_cell <- ccle_cell_lines[ccle_idx]
ccle_set <- cbind(ccle_exp[ccle_idx, ], auc = ccle_auc[ccle_idx, drug])
gdsc_idx <- which(!is.na(gdsc_auc[, drug]));        gdsc_cell <- gdsc_cell_lines$cell[gdsc_idx]
gdsc_set <- cbind(gdsc_exp[gdsc_idx, ], auc = gdsc_auc[gdsc_idx, drug])
# print(ccle_set[1:8, c(1:8, ncol(ccle_set))]);       print(gdsc_set[1:8, c(1:8, ncol(gdsc_set))])


## Feature selection...
# ccle_fs <- attrEval("auc", data = ccle_set, estimator = "RReliefFequalK")
# ccle_fs <- ccle_fs[order(ccle_fs, decreasing = TRUE)]
# relieff_fs$ccle$weights[names(ccle_fs), drug] <- ccle_fs;     relieff_fs$ccle$ranks[, drug] <- names(ccle_fs)
# 
# write.table(ccle_fs, file = sprintf("./Data/ccle_relieff_fs_%s.tsv", drug), sep = "\t", row.names = TRUE, col.names = FALSE)

# gdsc_fs <- attrEval("auc", data = gdsc_set, estimator = "RReliefFequalK")
# gdsc_fs <- gdsc_fs[order(gdsc_fs, decreasing = TRUE)]
# relieff_fs$gdsc$weights[names(gdsc_fs), drug] <- gdsc_fs;     relieff_fs$gdsc$ranks[, drug] <- names(gdsc_fs)
# 
# write.table(gdsc_fs, file = sprintf("./Data/gdsc_relieff_fs_%s.tsv", drug), sep = "\t", row.names = TRUE, col.names = FALSE)

ccle_fs <- relieff_fs$ccle$ranks[, drug];       #gdsc_fs <- relieff_fs$gdsc$ranks[, drug]


# K <- 100
# features <- names(gdsc_fs)[1:K]
# features <- names(ccle_fs)[1:K]
features <- ccle_fs[1:K]
# features <- gdsc_fs[1:K]


### Modeling...
# X1 <- ccle_set[, features];     Y1 <- ccle_set$auc          ## target = ccle, source = gdsc
# X2 <- gdsc_set[, features];     Y2 <- gdsc_set$auc
X1 <- gdsc_set[, features];     Y1 <- gdsc_set$auc          ## target = gdsc, source = ccle
X2 <- ccle_set[, features];     Y2 <- ccle_set$auc

# seed <- 7531

# DMTL.
# prediction <- DMTL(target_set = list("X" = X1, "y" = Y1), source_set = list("X" = X2, "y" = Y2), 
#                    method = "hist", size = 1e3, seed = seed, pred_all = TRUE) 
# Y1_pred <- prediction$mapped;     Y1_pred_src <- prediction$unmapped
prediction <- DMTL(target_set = list("X" = X1, "y" = Y1), source_set = list("X" = X2, "y" = Y2), 
                   use_density = FALSE, sample_size = 1e3, random_seed = seed, all_pred = TRUE) 
Y1_pred <- prediction$target;     Y1_pred_src <- prediction$source


# Baseline.
# source("RF_predict.R")
Y1_pred_base <- RF_predict(x_train = norm_data(X2), y_train = Y2, x_test = norm_data(X1), y_lims = c(0, 1), 
                           n_tree = 200, m_try = 0.4, random_seed = seed)


perf_mat <- as.data.frame(rbind(
  COR1   = c(DMTL = COR(Y1, Y1_pred, sw = 1),   DMTL_SS = COR(Y1, Y1_pred_src, sw = 1),   BL = COR(Y1, Y1_pred_base, sw = 1)),
  NRMSE1 = c(DMTL = NRMSE(Y1, Y1_pred, sw = 1), DMTL_SS = NRMSE(Y1, Y1_pred_src, sw = 1), BL = NRMSE(Y1, Y1_pred_base, sw = 1)),
  NMAE1  = c(DMTL = NMAE(Y1, Y1_pred, sw = 1),  DMTL_SS = NMAE(Y1, Y1_pred_src, sw = 1),  BL = NMAE(Y1, Y1_pred_base, sw = 1)), 
  COR2   = c(DMTL = COR(Y1, Y1_pred, sw = 2),   DMTL_SS = COR(Y1, Y1_pred_src, sw = 2),   BL = COR(Y1, Y1_pred_base, sw = 2)),
  NRMSE2 = c(DMTL = NRMSE(Y1, Y1_pred, sw = 2), DMTL_SS = NRMSE(Y1, Y1_pred_src, sw = 2), BL = NRMSE(Y1, Y1_pred_base, sw = 2)),
  NMAE2  = c(DMTL = NMAE(Y1, Y1_pred, sw = 2),  DMTL_SS = NMAE(Y1, Y1_pred_src, sw = 2),  BL = NMAE(Y1, Y1_pred_base, sw = 2))
))
if (verbose) {
print(perf_mat[1:3, ]);     #print(perf_mat[4:6, ])
}

# write.table(perf_mat, file = sprintf("./Data/results_ccle_prediction_%s.tsv", drug), sep = "\t",
#             row.names = TRUE, col.names = TRUE)
# write.table(perf_mat, file = sprintf("./Data/results_gdsc_prediction_%s.tsv", drug), sep = "\t",
#             row.names = TRUE, col.names = TRUE)

# results_ccle$DMTL[drug, ] <- perf_mat$DMTL;     results_ccle$DMTL_SS[drug, ] <- perf_mat$DMTL_SS
# results_ccle$BL[drug, ]   <- perf_mat$BL
results_gdsc$DMTL[drug, ] <- perf_mat$DMTL;     results_gdsc$DMTL_SS[drug, ] <- perf_mat$DMTL_SS
results_gdsc$BL[drug, ]   <- perf_mat$BL

# if (verbose) {
# #print(rbind(dim(X1), dim(X2)));     
# print(c(length(Y1), length(Y2)));      cat("\n")
# }

# }
# results_ccle <- as.data.frame(results_ccle)
results_gdsc <- as.data.frame(results_gdsc)
# }

# save(relieff_fs, file = "./Data/relieff_fs.RData")
# save(results_ccle, results_gdsc, file = "./Data/DMTL_results_ccle_gdsc.RData")
# write.table(relieff_fs$ccle$ranks, file = "./Data/ccle_drug_ranks_genes.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
# write.table(relieff_fs$ccle$weights, file = "./Data/ccle_drug_rrelieff_weights.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
# write.table(relieff_fs$gdsc$ranks, file = "./Data/gdsc_drug_ranks_genes.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
# write.table(relieff_fs$gdsc$weights, file = "./Data/gdsc_drug_rrelieff_weights.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


results_gdsc <- run(seed = 4444, K = 100, verbose = FALSE)

cc1 <- cbind(results_gdsc[, c("DMTL.CC1", "DMTL_SS.CC1", "BL.CC1")]);   cc1 <- rbind(cc1, mean = colMeans(cc1))
colnames(cc1) <- c("DMTL", "DMTL_SS", "BL");    #print(cc1)

cc2 <- cbind(results_gdsc[, c("DMTL.CC2", "DMTL_SS.CC2", "BL.CC2")]);   cc2 <- rbind(cc2, mean = colMeans(cc2))
colnames(cc2) <- c("DMTL", "DMTL_SS", "BL");    #print(cc2)

nrmse1 <- cbind(results_gdsc[, c("DMTL.NRMSE1", "DMTL_SS.NRMSE1", "BL.NRMSE1")])
nrmse1 <- rbind(nrmse1, mean = colMeans(nrmse1));     colnames(nrmse1) <- c("DMTL", "DMTL_SS", "BL");    #print(nrmse1)

nmae1 <- cbind(results_gdsc[, c("DMTL.NMAE1", "DMTL_SS.NMAE1", "BL.NMAE1")])
nmae1 <- rbind(nmae1, mean = colMeans(nmae1));        colnames(nmae1) <- c("DMTL", "DMTL_SS", "BL");    #print(nmae1)

print( rbind(cc1 = cc1["mean", ], cc2 = cc2["mean", ], nrmse1 = nrmse1["mean", ], nmae1 = nmae1["mean", ]) )





