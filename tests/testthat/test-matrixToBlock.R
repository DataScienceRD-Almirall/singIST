test_that("matrixToBlock processes superpathway.input correctly", {
    data("example_superpathway_input", package = "singIST")
    result <- matrixToBlock(example_superpathway_input)
    # Test that the function returns a list
    expect_type(result, "list")
    # Test that required elements are returned
    expect_true("block_predictor" %in% names(result))
    expect_true("matrix_response" %in% names(result))
    expect_true("block_dim" %in% names(result))
    expect_true("observed_gene_sets" %in% names(result))
    # Check if block_predictor is a matrix
    expect_class(result$block_predictor, "matrix")
    # Check if matrix_response is a matrix
    expect_class(result$matrix_response, "matrix")
    # Check if block_dim is a numeric vector
    expect_class(result$block_dim, "numeric")
    # Test for invalid object (non-superpathway.input object)
    expect_error(matrixToBlock("invalid_object"))
    # Test for case with empty gene sets
    # Create an object with empty gene sets
    data_with_empty_gene_sets <- example_superpathway_input
    data_with_empty_gene_sets@superpathway_info@gene_sets_celltype <- list()
    expect_error(matrixToBlock(data_with_empty_gene_sets))
})
