test_that("Check consistency of superpathway.input class slots", {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
    # Load pathway and superpathway.gene.sets objects
    my_pathway <- create_pathway(standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
                           dbsource = "KEGG",
                           collection = "c2",
                           subcollection = "CP")
    celltypes <- c("T-cell", "Dendritic Cell")
    my_superpathway <- create_superpathway(my_pathway, celltypes, list())
    # Load hyperparameters object
    quantile_comb_table <- base::as.matrix(
        RcppAlgos::permuteGeneral(seq(0.05, 0.95, by = 0.50)),
        ncol = length(celltypes)
        )
    outcome_type <- "binary"
    number_PLS <- as.integer(3)
    folds_CV <- as.integer(1)
    repetition_CV <- as.integer(1)
    my_hyperparameters <- create_hyperparameters(quantile_comb_table, outcome_type, number_PLS, folds_CV, repetition_CV)
    # Slots of superpathway input
    sample_id <- c("AD1", "AD2", "HC1", "HC2")
    sample_class <- c("AD", "AD", "HC", "HC")
    base_class <- "HC"
    target_class <- "AD"
    pseudobulk_lognorm <- matrix(rnorm(length(celltypes)*length(sample_id)),
    nrow = length(celltypes)*length(sample_id), ncol = length(celltypes))
    rownames(pseudobulk_lognorm) <- as.vector(t(outer(celltypes, sample_id,
    function(x, y) paste(x, y, sep = "_"))))

    # Expect pass
    expect_list(create_superpathway_input(
                     superpathway_info = my_superpathway,
                     hyperparameters_info = my_hyperparameters,
                     pseudobulk_lognorm = pseudobulk_lognorm,
                     sample_id = sample_id,
                     sample_class = sample_class,
                     base_class = base_class,
                     target_class = target_class))
})
