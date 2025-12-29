test_that("superpathway_recap works as expected", {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
    file <- system.file("extdata", "example_superpathway_fit_model.rda", package = "singIST")
    load(file)
    file <- system.file("extdata", "example_mapping_organism.rda", package = "singIST")
    load(file)
    model <- example_superpathway_fit_model
    mapped <- example_mapping_organism
    singIST_samples <- biological_link_function(mapped, model,
                                                exact = FALSE)$singIST_samples
    original <- derive_contributions(model, singIST_samples)
    derived <- derive_contributions(model,
                                    model$model_fit$predictor_block)
    result <- superpathway_recap(model, original$superpathway_score,
                                    derived$superpathway_score)
    # Checking if the result is a data frame
    expect_true(is.data.frame(result))
    # Checking if the required columns are present
    expect_true("pathway" %in% colnames(result))
    expect_true("recapitulation" %in% colnames(result))
})
