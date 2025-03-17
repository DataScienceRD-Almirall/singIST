test_that("permut_asmbplsda runs correctly", {
    data("example_superpathway_fit_model", package = "singIST")
    # Set up necessary parameters
    CV_error <- 0.05
    npermut <- 5
    Nc <- 1
    # Run the permutation test
    result <- permut_asmbplsda(example_superpathway_fit_model, npermut = npermut,
                               CV_error = CV_error, Nc = Nc)
    # Test that the result is a list and contains necessary elements
    expect_type(result, "list")
    expect_true("pvalue" %in% names(result))
    expect_true("IC" %in% names(result))
    expect_true("prct.Ychange.values" %in% names(result))
    # Check that p-value is a numeric value
    expect_type(result$pvalue, "double")
    # Check that the confidence intervals (IC) are also numeric
    expect_type(result$IC, "double")
    # Test for edge cases with different CV_error values (e.g., large CV_error)
    result_high_error <- permut_asmbplsda(example_superpathway_fit_model,
                                            npermut = npermut,
                                            CV_error = 2, Nc = Nc)
    expect_type(result_high_error$pvalue, "double")
})
