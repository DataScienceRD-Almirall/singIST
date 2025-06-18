test_that("biological_link_function runs end-to-end correctly", {
    file <- system.file("extdata", "example_mapping_organism.rda", package = "singIST")
    load(file)
    object <- example_mapping_organism
    file <- system.file("extdata", "example_superpathway_fit_model.rda", package = "singIST")
    load(file)
    model_object <- example_superpathway_fit_model
    result <- biological_link_function(object, model_object, exact = FALSE)
    # Check that the result is a list with expected elements
    expect_true(is.list(result))
    expect_true("orthologs" %in% names(result))
    expect_true("singIST_samples" %in% names(result))
    expect_true("FC" %in% names(result))
})
