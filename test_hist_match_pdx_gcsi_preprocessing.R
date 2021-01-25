# rm(list = ls())

## Set-up system path...
PATH <- if (Sys.getenv("USERNAME") == "SRDhruba") {
  "\\Users\\SRDhruba\\Dropbox (Personal)\\ResearchWork\\Rtest\\"
} else {
  sprintf("%s\\Dropbox\\ResearchWork\\Rtest\\", Sys.getenv("HOMEPATH"))
}
setwd(PATH);       cat("Current system path = ", getwd(), "\n")

DATA.PATH <- c(PDX  = sprintf("%s\\Dropbox (Texas Tech)\\Tumor_CL_TL_2018\\PDX_data_NatMed", Sys.getenv("HOMEPATH")),
               gCSI = sprintf("%s\\Dropbox (Texas Tech)\\Tumor_CL_TL_2018\\PDX_data_NatMed\\gCSI", Sys.getenv("HOMEPATH")))


## Packages...
library(openxlsx)


#### Functions...
printf <- function(..., end = "\n") {
  if ((nargs() > 1) & (grepl(list(...)[1], pattern = "%")))
    cat(sprintf(...), end)
  else
    cat(..., end)
}


## Read PDX data...
PDX.data <- list(
  RNA = read.xlsx(sprintf("%s\\PDX_gene_exp_log2rpkm_17470_genes_processed_SRD_07_06_2020.xlsx", DATA.PATH["PDX"]), sheet = 2),
  AUC = read.xlsx(sprintf("%s\\PDX_time_response_AUC_62_drugs_processed_SRD_02_27_2020.xlsx", DATA.PATH["PDX"]), sheet = 2)
)

gx.pdx  <- as.data.frame(t(PDX.data$RNA[, -1]));      colnames(gx.pdx) <- PDX.data$RNA[, 1]
auc.pdx <- as.data.frame(PDX.data$AUC[, -1], row.names = PDX.data$AUC[, 1])

all(rownames(gx.pdx) == rownames(auc.pdx))      # Check common samples


## Read gCSI data...
gCSI.data <- lapply(1:2, function(i) {
  read.xlsx(sprintf("%s\\gCSI_datasets_processed_for_modeling_24_Jun_2020.xlsx", DATA.PATH["gCSI"]), sheet = i)
})
names(gCSI.data) <- c("RNA", "AUC")

gx.cell  <- as.data.frame(t(gCSI.data$RNA[, -1]));     colnames(gx.cell) <- gCSI.data$RNA[, 1]
gx.cell  <- gx.cell[order(rownames(gx.cell)), order(colnames(gx.cell))]
auc.cell <- as.data.frame(gCSI.data$AUC[, -1], row.names = gCSI.data$AUC[, 1])
auc.cell <- auc.cell[order(rownames(auc.cell)), ]

all(rownames(gx.cell) == rownames(auc.cell))      # Check common samples


## Common datasets...
genes <- intersect(colnames(gx.pdx),  colnames(gx.cell));     printf("#common genes = %d", length(genes))
drugs <- intersect(tolower(colnames(auc.pdx)), tolower(colnames(auc.cell)));
printf("#common drugs = %d", length(drugs))

gx.pdx2 <- gx.pdx[, genes];     gx.cell2 <- gx.cell[, genes]














