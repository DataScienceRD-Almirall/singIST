test_that("superpathway_recap works as expected", {
    data("example_superpathway_fit_model", package = "singIST")
    data("example_mapping_organism", package = "singIST")
    model <- example_superpathway_fit_model
    mapped <- example_mapping_organism
    singIST_samples <- biological_link_function(mapped, model,
                                                exact = FALSE)$singIST_samples
    original <- derive_contributions(model, singIST_samples)
    derived <- derive_contributions(model,
                                    slot(model, "model_fit")$predictor_block)
    result <- superpathway_recap(model, original$superpathway_score,
                                    derived$superpathway_score)
    # Checking if the result is a data frame
    expect_true(is.data.frame(result))
    # Checking if the required columns are present
    expect_true("pathway" %in% colnames(result))
    expect_true("recapitulation" %in% colnames(result))
})
