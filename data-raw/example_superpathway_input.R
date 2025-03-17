## code to prepare `example_superpathway_input` object goes here
cytokine_pathway <- methods::new("pathway",
                                    standard_name = output_dataset$Pathway_name,
                                    dbsource = "KEGG",
                                    collection = "c2",
                                    subcollection = "CP")
# Gene set
gene_set <- c("CCL26", "TNFSF13", "HGF", "CCL3L1", "TNFSF12", "TNFRSF8",
            "TNFSF10", "CCL2", "TNFSF8", "TNFSF9", "CCL3", "TNFSF14", "IL21R",
            "CCL5", "CCL4", "CCL13", "IL17RB", "CCL11", "CCL8", "CCL7",
            "TNFRSF17", "TNFSF11", "IL23A", "CTF1", "FLT1", "CCL14", "FLT3",
            "CCL15", "FLT3LG", "FLT4", "IL22", "CD27", "CCL18", "CCL17",
            "CCL20", "CCL19", "CCL21", "CCL23", "CCL22", "TPO", "CCL16",
            "CSF2", "CSF1R", "TNFRSF13C", "IL13RA1", "CSF1", "PRLR", "PRL",
            "TNFRSF14", "IFNE", "CCL24", "CSF2RB", "CSF2RA", "TGFB2", "TGFB1",
            "CX3CL1", "IL23R", "XCL1", "CXCL5", "CXCL11", "CXCL6", "CCR2",
            "CCL25", "CSF3", "TGFBR2", "TGFBR1", "MET", "TGFB3", "CXCR3",
            "CSF3R", "CLCF1", "IL17B", "CXCL12", "CCR10", "OSMR", "XCR1",
            "IL20", "INHBE", "PLEKHO2", "CXCL10", "MPL", "CCL27", "INHBB",
            "INHBC", "IFNLR1", "INHBA", "KITLG", "TNFRSF19", "CXCR4", "IL1R2",
            "KIT", "CD70", "TNFRSF11B", "LEPR", "CXCR6", "RELT", "LEP",
            "IL19", "LTA", "IL17A", "IL18", "TNFRSF9", "OSM", "GDF5",
            "IL15RA", "IL15", "CXCL13", "LTB", "TNFRSF25", "LTBR", "IL22RA1",
            "CD40LG", "CD40", "PPBPP1", "GHR", "IL18RAP", "CXCL9",
            "IL17RA", "GH1", "LIF", "GH2", "LIFR", "IL18R1", "TNFSF4",
            "TNFRSF4", "CX3CR1", "CCL3L3", "PDGFRB", "EDA", "CXCL14", "IL1A",
            "EDAR", "IL24", "CNTF", "CNTFR", "TNFSF13B", "PDGFC", "VEGFA",
            "VEGFC", "VEGFB", "CCR9", "NGFR", "IL22RA2", "PDGFA", "ACVRL1",
            "PDGFB", "CCL4L2", "ACVR2B", "PDGFRA", "ACVR2A", "ACVR1B",
            "ACVR1", "VEGFD", "TNFRSF13B", "IL25", "TNFSF15", "EPO", "IL9R",
            "CXCL1", "TNF", "CXCR5", "IFNA5", "IFNA4", "IFNA2", "TNFRSF6B",
            "IFNA1", "CXCR2", "CXCR1", "IL9", "CCL28", "IL7R", "CXCL8",
            "IL12RB1", "IL12B", "IL13", "IL12RB2", "IL11RA", "IL12A", "KDR",
            "TNFRSF1B", "IFNA17", "TNFRSF1A", "IFNA21", "TNFRSF18", "IFNA6",
            "CCR3", "TNFRSF12A", "IFNA7", "CCR4", "IFNA8", "CCR5", "IFNA10",
            "CCR6", "IFNA13", "CCR7", "IFNA14", "CCR8", "IFNA16", "IL26",
            "CXCL16", "IL10", "EDA2R", "IL10RA", "IL10RB", "EPOR", "IL11",
            "IL3RA", "IL3", "IL2RG", "EGFR", "IL2RB", "CCR1", "XCL2",
            "IFNGR2", "EGF", "CCL1", "TNFRSF10A", "IFNG", "IFNGR1",
            "TNFRSF10D", "IFNB1", "TNFRSF11A", "TNFRSF10B", "IFNAR1",
            "TNFRSF10C", "IFNAR2", "IL1RAP", "IL1B", "AMHR2", "IL1R1", "AMH",
            "IL2RA", "IL20RB", "IL21", "TNFSF18", "IL2", "IL20RA", "BMPR2",
            "BMPR1A", "IL6R", "BMPR1B", "CRLF2", "IL6ST", "IFNK", "IL7",
            "PF4", "CXCL2", "CXCL3", "BMP2", "PF4V1", "BMP7", "IFNL3",
            "IFNL1", "TNFRSF21", "FAS", "IFNL2", "FASLG", "IFNW1", "PPBP",
            "TSLP", "IL4", "IL4R", "IL5", "IL5RA", "IL6"
)
cytokine_superpathway <- methods::new("superpathway.gene.sets",
                                        pathway_info = cytokine_pathway,
                                        celltypes = c("Dendritic Cells",
                                                        "Keratinocytes"),
                                      gene_sets_celltype = list(gene_set,
                                                                gene_set))
# Initialize hyperparameters
library(RcppAlgos)
quantile_comb_table <- as.matrix(
    RcppAlgos::permuteGeneral(seq(0.05, 0.95, by = 1),
                                m = length(cytokine_superpathway@celltypes),
                                TRUE))
my_hyperparameters <- methods::new("hyperparameters",
                                   quantile_comb_table = quantile_comb_table,
                                   outcome_type = "binary",
                                   number_PLS = as.integer(1),
                                   folds_CV = as.integer(1),
                                   repetition_CV = as.integer(1))
# Initialize superpathway input object
pseudobulk <- matrix(rnorm(12*length(gene_set)), nrow = 6*2, ncol = length(gene_set))
names1 <- expand.grid("Dendritic Cells", c(paste0("Target", 1:3), paste0("Base", 1:3)))
names2 <- expand.grid("Keratinocytes", c(paste0("Target", 1:3), paste0("Base", 1:3)))
rownames(pseudobulk) <- c(paste0(names1[,1], "_", names1[,2]), paste0(names2[,1], "_", names2[,2]))
colnames(pseudobulk) <- gene_set
cytokine_superpathway_input <- methods::new("superpathway.input",
                                      superpathway_info = cytokine_superpathway,
                                      hyperparameters_info = my_hyperparameters,
                                      pseudobulk_lognorm = pseudobulk,
                                      sample_id = c("Target1", "Target2", "Target3",
                                                    "Base1", "Base2", "Base3"),
                                      sample_class = c(rep("Target",3), rep("Base", 3)),
                                      base_class = "Base",
                                      target_class = "Target")
example_superpathway_input <- cytokine_superpathway_input
usethis::use_data(example_superpathway_input, overwrite = TRUE)
