test_that("Test multiple_fitOptimal function", {
    data("example_superpathway_input", package = "singIST")
    # Define example superpathway.input objects
    object1 <- object2 <- example_superpathway_input
    models <- list(object1, object2)
    # Test when using default arguments
    result <- multiple_fitOptimal(models)
    expect_equal(length(result), 2)
    expect_true(all(sapply(result, inherits, "superpathway.fit.model")))
    # Test when passing custom arguments
    result <- multiple_fitOptimal(models, type = c("jackknife", "subsampling"))
    expect_equal(length(result), 2)
    # Test when passing different parameters for each model
    result <- multiple_fitOptimal(models, nsubsampling = c(100, 50),
                                    npermut = c(10, 15))
    expect_equal(length(result), 2)
    # Test invalid case when number of arguments doesn't match number of models
    expect_error(multiple_fitOptimal(models, nsubsampling = c(100, 50, 25)))
})
