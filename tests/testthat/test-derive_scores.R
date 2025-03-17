test_that("derive_scores computes scores correctly", {
    data("example_superpathway_fit_model", package = "singIST")
    object <- example_superpathway_fit_model
    data <- biological_link_function(
        example_mapping_organism, example_superpathway_fit_model,
        exact = FALSE)$singIST_samples
    sample <- 2
    # Call the function
    scores_output <- derive_scores(object, data, sample)
    # Check if the output is a list
    expect_true(is.list(scores_output))
    # Check the structure of the output
    expect_named(scores_output, c("Delta", "delta", "gamma", "Gamma",
                                "Y_pred_num", "Y_weight"))
    # Check if Delta is a matrix and has expected dimensions
    expect_class(scores_output$Delta, "matrix")
    expect_equal(nrow(scores_output$Delta),
                length(object@superpathway_input@superpathway_info@celltypes))
    # Check Y_pred_num is a numeric vector
    expect_type(scores_output$Y_pred_num, "double")
    # Check if gamma is a matrix with expected dimensions
    expect_class(scores_output$gamma, "matrix")
    expect_equal(nrow(scores_output$gamma),
            length(object@superpathway_input@superpathway_info@celltypes))
})

test_that("derive_scores handles empty input correctly", {
    # Test with an empty matrix
    empty_data <- matrix(ncol = 0, nrow = 0)
    expect_error(derive_scores(example_superpathway_fit_model, empty_data, 1))
})
