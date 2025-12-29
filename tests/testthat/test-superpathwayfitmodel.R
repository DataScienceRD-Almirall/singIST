test_that("Check that slots of superpathway fit model class are consistent", {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
    file <- system.file("extdata", "example_superpathway_fit_model.rda", package = "singIST")
    load(file)
    expect_list(example_superpathway_fit_model)
})
