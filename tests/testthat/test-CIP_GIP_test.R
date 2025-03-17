test_that("CIP_GIP_test runs correctly", {
    data("example_superpathway_fit_model", package = "singIST")
    # Set up necessary parameters
    npermut <- 3
    maxiter <- 100
    # Run the CIP_GIP_test for jackknife method
    result_jackknife <- CIP_GIP_test(example_superpathway_fit_model,
                                        npermut = npermut,
                                        maxiter = maxiter, type = "jackknife")
    # Test that the result is a list with expected components
    expect_type(result_jackknife, "list")
    expect_true("variability_param" %in% names(result_jackknife))
    expect_true("NULL_CIP_GIP" %in% names(result_jackknife))
    expect_true("CIP_pvalue" %in% names(result_jackknife))
    expect_true("GIP_pvalue" %in% names(result_jackknife))
    # Check that CIP_pvalue and GIP_pvalue are lists of numeric values
    expect_type(result_jackknife$CIP_pvalue, "list")
    expect_type(result_jackknife$GIP_pvalue, "list")
    # Run the CIP_GIP_test for subsampling method
    result_subsampling <- CIP_GIP_test(example_superpathway_fit_model,
                                        npermut = npermut,
                                        maxiter = maxiter, type = "subsampling")
    # Test that the result for subsampling is similar to jackknife
    expect_type(result_subsampling, "list")
    expect_true("variability_param" %in% names(result_subsampling))
    expect_true("NULL_CIP_GIP" %in% names(result_subsampling))
    expect_true("CIP_pvalue" %in% names(result_subsampling))
    expect_true("GIP_pvalue" %in% names(result_subsampling))
    # Check that p-values are properly calculated (numeric and not NA)
    expect_type(result_subsampling$CIP_pvalue, "list")
    expect_type(result_subsampling$GIP_pvalue, "list")
    expect_true(all(!is.na(unlist(result_subsampling$CIP_pvalue))))
    expect_true(all(!is.na(result_subsampling$GIP_pvalue)))
    # Test invalid type input
    expect_error(CIP_GIP_test(example_superpathway_fit_model, npermut = npermut,
                              maxiter = maxiter, type = "invalid_type"))
})

