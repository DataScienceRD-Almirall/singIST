test_that("Test multiple_fitOptimal function", {
    file <- system.file("extdata", "example_superpathway_input.rda", package = "singIST")
    load(file)
    # Define example superpathway.input objects
    object1 <- object2 <- example_superpathway_input
    models <- list(object1, object2)
    # Test when using default arguments
    result <- multiple_fitOptimal(models, exact = FALSE)
    expect_equal(length(result), 2)
    expect_true(all(sapply(result, inherits, "superpathway.fit.model")))
    # Test when passing custom arguments
    result <- multiple_fitOptimal(models, type = c("jackknife", "subsampling"),
                                  exact = FALSE)
    expect_equal(length(result), 2)
    # Test invalid case when number of arguments doesn't match number of models
    expect_error(multiple_fitOptimal(models, nsubsampling = c(100, 50, 25), exact = FALSE))
})
