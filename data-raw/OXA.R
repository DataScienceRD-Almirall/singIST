## code to prepare `OXA_IMQ` dataset goes here
library(readr)
library(Seurat)
# Load Seurat object with the disease model data
raw_OXA_IMQ <- readRDS("data-raw/OXA_IMQ.rds")
raw_OXA_IMQ <- Seurat::UpdateSeuratObject(raw_OXA_IMQ)
# Filter for Oxazolone model and cell types to model
OXA <- raw_OXA_IMQ[, raw_OXA_IMQ$stim %in% c("OXA", "ETOH")]
OXA <- OXA[, OXA$identities %in% c("cDC2", "cDC1", "migratory DCs",
                                    "Keratinocytes", "LC", "DETC", "dγdT", "T")]
# Subsample cells and genes of the object to reduce memory 
set.seed(123)
n_cells <- round(ncol(OXA)/22, 0)
n_genes <- round(nrow(OXA)/3, 0)
subsample_cells <- sample(colnames(OXA), n_cells)
subsample_genes <- sample(rownames(OXA), n_genes)
OXA <- subset(OXA, features = subsample_genes)
OXA <- subset(OXA, cells = subsample_cells)
usethis::use_data(OXA, compress = "xz", overwrite = TRUE)
