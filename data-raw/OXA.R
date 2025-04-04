## code to prepare `OXA_IMQ` dataset goes here
library(readr)
library(Seurat)
librar(dplyr)
# Load Seurat object with the disease model data
raw_OXA_IMQ <- readRDS("data-raw/OXA_IMQ.rds")
raw_OXA_IMQ <- Seurat::UpdateSeuratObject(raw_OXA_IMQ)
# Filter for Oxazolone model and cell types to model
OXA <- raw_OXA_IMQ[, raw_OXA_IMQ$stim %in% c("OXA", "ETOH")]
OXA <- OXA[, OXA$identities %in% c("cDC2", "cDC1", "migratory DCs",
                                    "Keratinocytes", "LC", "DETC", "dγdT", "T")]
# Subsample cells and genes of the object to reduce memory 
# Sample 20 cells per cell type + condition combination
meta <- OXA@meta.data  %>% mutate(cell_id = rownames(.))
sampled_cells <- meta %>% group_by(identities, stim)  %>% slice_sample(n = 20) %>% pull(cell_id)
sampled_genes <- sample(rownames(OXA), round(nrow(OXA)/2, 0))
OXA <- subset(OXA, features = sampled_genes)
OXA <- subset(OXA, cells = sampled_cells)
usethis::use_data(OXA, compress = "xz", overwrite = TRUE)
