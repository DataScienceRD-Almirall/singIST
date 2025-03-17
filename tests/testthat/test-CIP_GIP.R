test_that("CIP_GIP computes correctly", {
    data("example_superpathway_fit_model", package = "singIST")
    result <- CIP_GIP(example_superpathway_fit_model)
    # Test that result is a list with GIP and CIP components
    expect_type(result, "list")
    expect_true("GIP" %in% names(result))
    expect_true("CIP" %in% names(result))
    # Check that CIP is a matrix with cell types as row names
    expect_class(result$CIP, "matrix")
    # Test invalid object
    expect_error(CIP_GIP("invalid_object"))
    # Check that CIP and GIP sum to 1
    expect_equal(sum(result$CIP), 1)
    expect_equal(
        sum(result$GIP[[sample(length(result$GIP), size = 1)]][,1]), 1 )
})
