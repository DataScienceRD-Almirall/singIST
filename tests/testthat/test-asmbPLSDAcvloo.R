test_that("asmbPLSDA.cv.loo runs correctly", {
    data("example_superpathway_input", package = "singIST")
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
                                cores = NULL, expected.measure.increase = 0.005,
                                center = TRUE, scale = TRUE, maxiter = 100)
    # Test that result is a list and contains required components
    expect_type(result, "list")
    expect_true("quantile_table_CV" %in% names(result))
    expect_true("optimal_nPLS" %in% names(result))
    # Check the quantile_table_CV is a matrix
    expect_class(result$quantile_table_CV, "matrix")
    # Check optimal_nPLS is a numeric value
    expect_type(result$optimal_nPLS, "double")
    # Test parallel execution (assuming you have cores available)
    result_parallel <- asmbPLSDA.cv.loo(
        X.matrix, Y.matrix, PLS_term = 1, X.dim = X.dim,
        quantile.comb.table = quantile.comb.table, outcome.type = outcome.type,
        Method = NULL, measure = "B_accuracy", parallel = TRUE, cores = 2,
        expected.measure.increase = 0.005, center = TRUE, scale = TRUE,
        maxiter = 100)
    # Ensure parallel version returns similar results
    expect_equal(result$quantile_table_CV, result_parallel$quantile_table_CV)
})
