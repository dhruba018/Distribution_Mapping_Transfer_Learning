setwd("C:/Users/SRDhruba/Dropbox (Personal)/ResearchWork/Rtest/")

library(openxlsx)


### Read data...
ccle_exp_data <- read.xlsx("./Data/CCLE_gene_expression_Oct_13_SRD.xlsx", sheet = 1, check.names = FALSE, 
                           rowNames = TRUE, colNames = TRUE)
ccle_auc_data <- read.xlsx("./Data/CCLE_AUC_24_Drugs_processed_SRD_01_19_2021.xlsx", sheet = 1, check.names = FALSE, 
                           rowNames = TRUE, colNames = TRUE)
gdsc_exp_data <- read.xlsx("./Data/gdsc_v7_gene_expression_processed_SRD_09_27_2018.xlsx", sheet = 1, check.names = FALSE, 
                           rowNames = TRUE, colNames = TRUE)
gdsc_auc_data <- read.xlsx("./Data/gdsc_v7_dose_response_processed_SRD_01_19_2021.xlsx", sheet = "AUC", check.names = FALSE, 
                           rowNames = TRUE, colNames = TRUE)
gdsc_drug_info1 <- read.xlsx("./Data/gdsc_v7_dose_response_processed_SRD_01_19_2021.xlsx", sheet = "Drug_info", 
                             check.names = FALSE, rowNames = FALSE, colNames = TRUE)
gdsc_drug_info2 <- read.xlsx("./Data/gdsc_v7_Screened_Compounds.xlsx", sheet = 1, check.names = FALSE, 
                            rowNames = FALSE, colNames = TRUE)


### Process data...
check_dim <- function(x) rbind(all = c(row = nrow(x), col = ncol(x)), unique = c(row = length(unique(rownames(x))), 
                                                                                 col = length(unique(colnames(x)))))

## CCLE...
# Expression.
ccle_exp <- t(ccle_exp_data);     print(check_dim(ccle_exp))
rows_rep <- rownames(ccle_exp)[duplicated(rownames(ccle_exp))]
for (i in rows_rep) {
  ccle_exp[i, ] <- colMeans(ccle_exp[which(rownames(ccle_exp) == i), ])
}
ccle_exp <- ccle_exp[which(!duplicated(rownames(ccle_exp))), ]
print(check_dim(ccle_exp))
dimnames(ccle_exp) <- lapply(dimnames(ccle_exp), tolower)
colnames(ccle_exp) <- gsub(colnames(ccle_exp), pattern = "-", replacement = "_")
colnames(ccle_exp) <- gsub(colnames(ccle_exp), pattern = ".", replacement = "_", fixed = TRUE)
# print(ccle_exp[1:8, 1:8])


# AUC response.
ccle_auc <- as.matrix(ccle_auc_data);     print(check_dim(ccle_auc))
dimnames(ccle_auc) <- lapply(dimnames(ccle_auc), tolower)
colnames(ccle_auc) <- gsub(colnames(ccle_auc), pattern = "-", replacement = "_")
colnames(ccle_auc) <- gsub(colnames(ccle_auc), pattern = ".", replacement = "_", fixed = TRUE)
# print(ccle_auc[1:8, 1:8])


# Common cell line data.
ccle_cell_lines <- sort(intersect(rownames(ccle_exp), rownames(ccle_auc)))
# print(ccle_cell_lines[1:8])
ccle_exp <- ccle_exp[ccle_cell_lines, ];      ccle_auc <- ccle_auc[ccle_cell_lines, ]
# print(ccle_exp[1:8, 1:8]);                    print(ccle_auc[1:8, 1:8])


## GDSC...
# Expression.
gdsc_exp <- t(gdsc_exp_data);     print(check_dim(gdsc_exp))
gdsc_exp <- gdsc_exp[-which(rownames(gdsc_exp) == "TT_aero_dig_tract"), ]       # Remove data for a cell line with wrong tissue
print(check_dim(gdsc_exp))
rows_rep <- rownames(gdsc_exp)[which(duplicated(rownames(gdsc_exp)))]
for (i in rows_rep) {
  gdsc_exp[i, ] <- colMeans(gdsc_exp[which(rownames(gdsc_exp) == i), ])
}
gdsc_exp <- gdsc_exp[which(!duplicated(rownames(gdsc_exp))), ]
print(check_dim(gdsc_exp))
dimnames(gdsc_exp) <- lapply(dimnames(gdsc_exp), tolower)
colnames(gdsc_exp) <- gsub(colnames(gdsc_exp), pattern = "-", replacement = "_")
colnames(gdsc_exp) <- gsub(colnames(gdsc_exp), pattern = ".", replacement = "_", fixed = TRUE)
# print(gdsc_exp[1:8, 1:8])


# AUC response.
gdsc_auc <- as.matrix(gdsc_auc_data);     print(check_dim(gdsc_auc))
dimnames(gdsc_auc) <- lapply(dimnames(gdsc_auc), tolower)
colnames(gdsc_auc) <- gsub(colnames(gdsc_auc), pattern = "-", replacement = "_")
colnames(gdsc_auc) <- gsub(colnames(gdsc_auc), pattern = ".", replacement = "_", fixed = TRUE)
# print(gdsc_auc[1:8, 1:8])


# Common cell line data.
gdsc_cell_lines <- lapply(rownames(gdsc_exp), function(x) {                     # Get cell line & tissue info
  x <- strsplit(x, split = "_", fixed = TRUE)[[1]]
  x <- c(x[1], if (length(x) > 2) Reduce(x[-1], f = function(...) paste(..., sep = "_")) else x[2])
  x  
})
gdsc_cell_lines <- as.data.frame(matrix(unlist(gdsc_cell_lines), ncol = 2, byrow = TRUE), stringsAsFactors = FALSE)
dimnames(gdsc_cell_lines) <- list(gdsc_cell_lines[, 1], c("cell", "tissue"))
rownames(gdsc_exp) <- gdsc_cell_lines$cell                                      # Replace expression matrix rows with cell names 
# print(gdsc_cell_lines[1:8, ]);                    print(gdsc_exp[1:8, 1:8])

gdsc_cell_lines <- gdsc_cell_lines[sort(intersect(gdsc_cell_lines$cell, rownames(gdsc_auc))), ]
# print(gdsc_cell_lines[1:8, ])
gdsc_exp <- gdsc_exp[gdsc_cell_lines$cell, ];     gdsc_auc <- gdsc_auc[gdsc_cell_lines$cell, ]
# print(gdsc_exp[1:8, 1:8]);                        print(gdsc_auc[1:8, 1:8])



### Get common datasets...
## Expression.
genes       <- sort(intersect(colnames(ccle_exp), colnames(gdsc_exp)))
cell_lines  <- sort(intersect(ccle_cell_lines, gdsc_cell_lines$cell))
tissue_info <- gdsc_cell_lines[cell_lines, ]
# print(genes[1:8]);            print(length(genes))
# print(cell_lines[1:8]);       print(length(cell_lines))
# print(tissue_info[1:8, ]);    print(dim(tissue_info))

ccle_exp <- ccle_exp[, genes];     gdsc_exp <- gdsc_exp[, genes]
# ccle_exp <- ccle_exp[cell_lines, genes];     gdsc_exp <- gdsc_exp[cell_lines, genes]
# print(ccle_exp[1:8, 1:8]);                   print(gdsc_exp[1:8, 1:8])


## AUC response.
# Get drug info.
info_idx  <- match(gdsc_drug_info1$DRUG_ID, table = gdsc_drug_info2$DRUG_ID)
drug_info <- cbind(gdsc_drug_info1, synonyms = gdsc_drug_info2$SYNONYMS[info_idx])
colnames(drug_info) <- tolower(colnames(drug_info))
drug_info <- drug_info[, c("drug_name", "synonyms", "putative_target")]
drug_info$drug_name <- tolower(drug_info$drug_name);      drug_info$synonyms <- tolower(drug_info$synonyms)
drug_info$drug_name <- gsub(drug_info$drug_name, pattern = "-", replacement = "_")
drug_info$drug_name <- gsub(drug_info$drug_name, pattern = ".", replacement = "_", fixed = TRUE)
drug_info$drug_name <- gsub(drug_info$drug_name, pattern = " ", replacement = "_")
drug_info$synonyms  <- gsub(drug_info$synonyms, pattern = "-", replacement = "_")
drug_info$synonyms  <- gsub(drug_info$synonyms, pattern = ".", replacement = "_", fixed = TRUE)
drug_info$synonyms  <- sapply(drug_info$synonyms, function(x) {
  x <- strsplit(x, split = ", ", fixed = TRUE)[[1]]
  x <- sapply(x, function(xx) gsub(xx, pattern = " ", replacement = "_"))
  x <- Reduce(x, f = function(...) paste(..., sep = ", "))
})
# print(drug_info[1:8, ])

# Find a drug by name or synonym.
get_drug_idx <- function(drugs) {
  indices <- sapply(drugs, function(dd) {
    idx <- which(drug_info$drug_name == dd)
    idx <- if (!length(idx)) grep(dd, drug_info$drug_name) else idx
    idx <- if (!length(idx)) grep(dd, drug_info$synonyms) else idx
    idx
  })
  indices <- as.numeric(indices);     names(indices) <- drugs
  indices <- indices[!is.na(indices)]
  indices
}

drug_info <- drug_info[get_drug_idx(colnames(gdsc_auc)), ];     rownames(drug_info) <- colnames(gdsc_auc)
# print(drug_info[1:8, ])


# Common drugs.
drug_idx  <- get_drug_idx(drugs = colnames(ccle_auc));            drugs <- names(drug_idx)
drug_info <- cbind(drug_info[drug_idx, ], common_name = drugs);   rownames(drug_info) <- drugs

ccle_auc <- ccle_auc[, drugs];        gdsc_auc <- gdsc_auc[, drug_idx];         colnames(gdsc_auc) <- drugs
# print(ccle_auc[1:8, 1:8]);            print(gdsc_auc[1:8, 1:8])


### Final files => ccle_exp, gdsc_exp, ccle_auc, gdsc_auc, genes, ccle_cell_lines, gdsc_cell_lines, 
### cell_lines, tissue_info, drugs, drug_info


### Save files...
save(ccle_exp, gdsc_exp, ccle_auc, gdsc_auc, ccle_cell_lines, gdsc_cell_lines, cell_lines, tissue_info, genes, drugs, drug_info, 
     file = "./Data/ccle_gdsc_common_processed_SRD_19_Jan_2021.RData")
write.table(ccle_exp, file = "./Data/ccle_exp_common_processed_19_Jan_2021.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(ccle_auc, file = "./Data/ccle_auc_common_processed_19_Jan_2021.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(gdsc_exp, file = "./Data/gdsc_exp_common_processed_19_Jan_2021.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(gdsc_auc, file = "./Data/gdsc_auc_common_processed_19_Jan_2021.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(tissue_info, file = "./Data/ccle_gdsc_cell_lines_tissue_common_processed_19_Jan_2021.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE)
write.table(drug_info, file = "./Data/ccle_gdsc_drug_info_common_processed_19_Jan_2021.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE)







