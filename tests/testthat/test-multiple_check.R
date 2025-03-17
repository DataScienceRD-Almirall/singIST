test_that("Test multiple_check function", {
    # Test when parameter is NULL
    parameter <- NULL
    objectLength <- 3
    result <- multiple_check(parameter, objectLength)
    expect_equal(length(result), objectLength)
    expect_true(all(sapply(result, is.null)))
    # Test when parameter is a single value
    parameter <- "test"
    objectLength <- 5
    result <- multiple_check(parameter, objectLength)
    expect_equal(length(result), objectLength)
    expect_true(all(sapply(result, function(x) x == "test")))
    # Test when parameter length matches object length
    parameter <- c(1, 2, 3, 4, 5)
    objectLength <- 5
    result <- multiple_check(parameter, objectLength)
    expect_equal(length(result), objectLength)
    expect_equal(result, list(1, 2, 3, 4, 5))
    # Test when parameter length does not match object length
    parameter <- c(1, 2)
    objectLength <- 5
    expect_error(multiple_check(parameter, objectLength))
})
