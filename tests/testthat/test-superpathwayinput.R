test_that("Check consistency of superpathway.input class slots", {
    # Load pathway and superpathway.gene.sets objects
    my_pathway <- new("pathway",
        standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
        dbsource = "KEGG",
        collection = "c2",
        subcollection = "CP")
    celltypes <- c("T-cell", "Dendritic Cell")
    my_superpathway <- new("superpathway.gene.sets",
                            pathway_info = my_pathway,
                            celltypes = celltypes)
    # Load hyperparameters object
    quantile_comb_table <- base::as.matrix(
        RcppAlgos::permuteGeneral(seq(0.05, 0.95, by = 0.50)),
        ncol = length(celltypes)
        )
    outcome_type <- "binary"
    number_PLS <- as.integer(3)
    folds_CV <- as.integer(1)
    repetition_CV <- as.integer(1)
    my_hyperparameters <- new("hyperparameters",
                                quantile_comb_table = quantile_comb_table,
                                outcome_type = outcome_type,
                                number_PLS = number_PLS,
                                folds_CV = folds_CV
                                )
    # Slots of superpathway input
    sample_id <- c("AD1", "AD2", "HC1", "HC2")
    sample_class <- c("AD", "AD", "HC", "HC")
    base_class <- "HC"
    target_class <- "AD"
    pseudobulk_lognorm <- matrix(rnorm(length(celltypes)*length(sample_id)),
    nrow = length(celltypes)*length(sample_id), ncol = length(celltypes))
    rownames(pseudobulk_lognorm) <- as.vector(t(outer(celltypes, sample_id,
    function(x, y) paste(x, y, sep = "_"))))

    # Expect error due to incorrect rownames not matching cell type and
    # sample id specified
    pseudobulk_lognorm_error <- pseudobulk_lognorm
    rownames(pseudobulk_lognorm_error)[4] <- "T-cell"
    expect_error(new("superpathway.input",
                        superpathway_info = my_superpathway,
                        hyperparameters_info = my_hyperparameters,
                        pseudobulk_lognorm = pseudobulk_lognorm_error,
                        sample_id = sample_id,
                        sample_class = sample_class,
                        base_class = base_class,
                        target_class = target_class))
    # Expect pass
    expect_class(new("superpathway.input",
                     superpathway_info = my_superpathway,
                     hyperparameters_info = my_hyperparameters,
                     pseudobulk_lognorm = pseudobulk_lognorm,
                     sample_id = sample_id,
                     sample_class = sample_class,
                     base_class = base_class,
                     target_class = target_class), "superpathway.input")
})
