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
usethis::use_data(OXA, compress = "xz", overwrite = TRUE)
