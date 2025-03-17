test_that("wilcox_CIP_GIP computes p-value correctly", {
    # Generate example distributions
    ref_distr <- rnorm(100, mean = 30, sd = 2)
    null_distr <- rnorm(100, mean = 0, sd = 1)
    result <- wilcox_CIP_GIP(ref_distr, null_distr)
    # Test that the result is a numeric p-value
    expect_type(result, "double")
    expect_length(result, 1)
    # Test edge case with greater distributions (expected p-value should be 1)
    ref_distr_identical <- rnorm(10000, mean = 1, sd = 1)
    null_distr_identical <- rnorm(10000, mean = 100, sd = 1)
    result_identical <- wilcox_CIP_GIP(ref_distr_identical, null_distr_identical)
    expect_equal(result_identical, 1, tolerance = 1e-3)
    # Test with different distributions (expected p-value should be < 0.05)
    result_diff <- wilcox_CIP_GIP(ref_distr, null_distr)
    expect_true(result_diff < 0.05)
})
