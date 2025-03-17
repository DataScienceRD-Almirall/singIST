test_that("fitOptimal.superpathway.input runs correctly", {
    data("example_superpathway_input", package = "singIST")
    # Set up necessary parameters
    npermut <- 5
    type <- "jackknife"
    nsubsampling <- 10
    # Run the optimal fit for the model
    result <- fitOptimal.superpathway.input(example_superpathway_input,
                                            npermut = npermut, type = type,
                                            nsubsampling = nsubsampling)
    # Test that the result is a 'superpathway.fit.model' object
    expect_class(result, "superpathway.fit.model")
    # Check that the hyperparameters are correctly returned
    expect_class(result@hyperparameters_fit, "hyperparameters")
    expect_type(result@hyperparameters_fit@number_PLS, "integer")
    expect_class(result@hyperparameters_fit@quantile_comb_table, "matrix")
    expect_type(result@hyperparameters_fit@folds_CV, "integer")
    # Test that model fit contains necessary information
    expect_true("predictor_block" %in% names(result@model_fit))
    expect_true("response_matrix" %in% names(result@model_fit))
    expect_true("observed_gene_sets" %in% names(result@model_fit))
    expect_true("asmbPLS-DA" %in% names(result@model_fit))
    # Check that the model fit results contain the expected PLS component data
    expect_type(result@model_fit$`asmbPLS-DA`, "list")
    # Run the same function with different parameters
    result_high_permut <- fitOptimal.superpathway.input(
        example_superpathway_input, npermut = 50, type = "subsampling",
        nsubsampling = 20)
    # Test that the result with higher permutations is still valid
    expect_class(result_high_permut, "superpathway.fit.model")
    expect_type(result_high_permut@hyperparameters_fit@number_PLS, "integer")
    expect_true("asmbPLS-DA" %in% names(result_high_permut@model_fit))
    # Check that global significance full and CIP.GIP significance full work
    # correctly
    result_with_full_significance <- fitOptimal.superpathway.input(
        example_superpathway_input, npermut = 5, type = "jackknife",
        global_significance_full = TRUE, CIP.GIP_significance_full = TRUE)
    # Validate that full significance results are returned in the object
    expect_class(result_with_full_significance, "superpathway.fit.model")
})

