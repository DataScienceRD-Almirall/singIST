## code to prepare `OXA_IMQ` dataset goes here
library(readr)
library(Seurat)
library(dplyr)
# Load Seurat object with the disease model data
file <- system.file("extdata", "example_mapping_organism.rda", package = "singIST")
load(file)
raw_OXA_IMQ <- readRDS("data-raw/OXA_IMQ.rds")
raw_OXA_IMQ <- Seurat::UpdateSeuratObject(raw_OXA_IMQ)
# Filter for Oxazolone model and cell types to model
OXA <- raw_OXA_IMQ[, raw_OXA_IMQ$stim %in% c("OXA", "ETOH")]
OXA <- OXA[, OXA$identities %in% c("cDC2", "cDC1", "migratory DCs",
                                    "Keratinocytes", "LC", "DETC", "T")]
# genes to retain
retain <- 
    c("IL3", "IL13", "IFNG", "ANPEP", "CSF2", "IL5", "CD7", "IFNA1",
      "CD5", "IL4", "CD33", "IFNB1", "IL10", "ITGAX", "CD40", "CD2",
      "TLR7", "IL3", "IL13", "IFNG", "ANPEP", "CSF2", "IL5", "CD7",
      "IFNA1", "CD5", "IL4", "CD33", "IFNB1", "IL10", "ITGAX", "CD40",
      "CD2", "TLR7", "IL3", "IL13", "IFNG", "ANPEP", "CSF2", "IL5",
      "CD7", "IFNA1", "CD5", "IL4", "CD33", "IFNB1", "IL10", "ITGAX",
      "CD40", "CD2", "TLR7", "IL2RG", "TNFSF10", "IL1RAP", "VEGFC", "ACVR2B",
      "EGF", "IL7", "CCR8", "CXCL14", "CCL15", "CXCR5", "PRL", "TNFRSF10A",
      "LIF", "IL12A", "CXCR4", "LTB", "TGFB3", "HGF", "BMPR1B", "INHBC",
      "TNFSF8", "CCL21", "IFNA5", "CCL1", "BMP7", "IL1R2", "IL17B", "IL3",
      "CSF3", "CXCL3", "IL17A", "IL17RB", "PDGFC", "IL22", "IL10", "ACVR2A",
      "FASLG", "TNFRSF12A", "IFNA17", "TNFRSF4", "GDF5", "TNFSF9", "CXCL10", "CTF1",
      "IFNE", "IFNA10", "PDGFRB", "PPBP", "CCR7", "IL18", "TNFSF10", "IFNL1",
      "CXCR6", "FLT4", "CXCL16", "CCL11", "IL11", "CCL15", "FAS", "AMH",
      "CCL19", "CXCL12", "TNFSF18", "IL12B", "CD70", "PF4", "IL2RA", "ACVR2A",
      "CCL24", "FASLG", "IFNA14", "CD27", "IFNA1", "KDR", "IFNA5", "TNFRSF10A",
      "IL23R", "CRLF2", "TNFSF12", "CCL4L2", "TGFBR2", "IL10RA", "IFNA7", "IFNK",
      "TNFRSF13B", "PRLR", "IL21", "TNFRSF6B", "CXCL11", "TNFRSF9", "TNFRSF21", "IL1R1",
      "PDGFRA", "VEGFD", "IFNA16", "IL3RA", "TGFBR1", "TNF", "IL18RAP", "CXCL8",
      "IL20", "IL11RA", "TNFSF14", "BMP7", "CCL1", "TNFRSF1B", "IL12A", "IL10",
      "IFNL3", "IL20RA", "TNFRSF19", "TNFRSF14", "CXCL1", "IFNA4", "CX3CR1", "IFNAR1",
      "TNFRSF4", "CCL5", "IL17B", "IL21R", "TNFRSF18", "EGF", "CCL21", "BMPR2",
      "IL7", "CCL25", "FLT3LG", "TNFRSF1A", "CNTF", "ACVR2B", "KIT", "IFNGR1",
      "CSF1R", "BMPR1B", "IL22", "OSMR", "GDF5", "INHBC", "IL26", "PF4V1",
      "IFNA8", "IL17RA", "CXCR3", "IL15RA", "TNFSF9", "TGFB1", "CXCL6", "AMH",
      "IFNLR1", "FLT3", "CCL15", "PRLR", "IFNA14", "CXCR4", "IL2RB", "CXCL16",
      "IL10", "CCL13", "IL21R", "CSF2RA", "EPO", "IFNA10", "XCL2", "IL4",
      "EGFR", "VEGFD", "CXCL6", "IL17A", "CCL4L2", "TSLP", "NGFR", "IL13",
      "CCL18", "IFNB1", "CLCF1", "LEPR", "CCL22", "IFNL3", "IL11", "TNF",
      "FLT1", "TNFRSF17", "IL23R", "PLEKHO2", "KDR", "IL19", "IFNGR1", "IL2RB",
      "CCR8", "IL1A", "TNFRSF10C", "GHR", "IFNW1", "CXCL8", "IFNA6", "EGF",
      "INHBC", "CCR5", "ACVR2B", "GDF5", "VEGFA", "TNFSF14", "IL22RA1", "INHBA",
      "IL10RA", "CSF1", "BMPR1A", "CXCR3", "BMPR1B", "IFNLR1", "IL4R", "CCL21",
      "LIF", "CXCL11", "BMP2", "GH2", "IL6R", "CCL27", "TNFSF14", "IFNB1",
      "CCR6", "CNTFR", "CXCR6", "IL22RA2", "EDA2R", "LTB", "IL7R", "CSF3",
      "IFNAR1", "VEGFC", "INHBC", "CCL28", "FAS", "CCL20", "AMH", "CCR1",
      "VEGFB", "CCL21", "CCR3", "CLCF1", "IL6", "TNFRSF10A", "BMP7", "LEPR")
genes_formateados <- stringr::str_to_title(tolower(retain))
# Subsample cells and genes of the object to reduce memory 
# Sample 20 cells per cell type + condition combination
meta <- OXA@meta.data  %>% mutate(cell_id = rownames(.))
#sampled_cells <- meta %>% group_by(identities, stim)  %>% slice_sample(n = 20) %>% pull(cell_id)
#sampled_genes <- sample(rownames(OXA), round(nrow(OXA)/2, 0))
OXA <- subset(OXA, features = genes_formateados)
#OXA <- subset(OXA, cells = sampled_cells)
OXA <- Seurat::DietSeurat(OXA, layer = c("counts", "data"), dimreducs = NULL, graphs = NULL)
usethis::use_data(OXA, compress = "xz", overwrite = TRUE)
