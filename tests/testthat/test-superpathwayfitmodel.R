test_that("Check that slots of superpathway fit model class are consistent", {
    file <- system.file("extdata", "example_superpathway_fit_model.rda", package = "singIST")
    load(file)
    expect_class(example_superpathway_fit_model, "superpathway.fit.model")
    expect_error(
        example_superpathway_fit_model@hyperparameters_fit@folds_CV <- -1)
})
