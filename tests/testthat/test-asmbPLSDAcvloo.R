test_that("asmbPLSDA.cv.loo runs correctly", {
    file <- system.file("extdata", "example_superpathway_input.rda", package = "singIST")
    load(file)
    matrices <- matrixToBlock(example_superpathway_input)
    X.matrix <- matrices$block_predictor
    Y.matrix <- matrices$matrix_response
    X.dim <- matrices$block_dim
    quantile.comb.table <- slot(slot(example_superpathway_input,
                                "hyperparameters_info"),
                                "quantile_comb_table")
    outcome.type <- slot(slot(example_superpathway_input,
                            "hyperparameters_info"),
                            "outcome_type")
    result <- asmbPLSDA.cv.loo(X.matrix, Y.matrix, PLS_term = 1, X.dim = X.dim,
                                quantile.comb.table = quantile.comb.table,
                                outcome.type = outcome.type, Method = NULL,
                                measure = "B_accuracy", parallel = FALSE,
                                expected.measure.increase = 0.005,
                                center = TRUE, scale = TRUE, maxiter = 100)
    # Test that result is a list and contains required components
    expect_type(result, "list")
    expect_true("quantile_table_CV" %in% names(result))
    expect_true("optimal_nPLS" %in% names(result))
    # Check the quantile_table_CV is a matrix
    expect_class(result$quantile_table_CV, "matrix")
    # Check optimal_nPLS is a numeric value
    expect_type(result$optimal_nPLS, "double")
})
