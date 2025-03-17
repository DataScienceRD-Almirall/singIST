test_that("biological_link_function runs end-to-end correctly", {
    data("example_mapping_organism", package = "singIST")
    object <- example_mapping_organism
    model_object <- example_superpathway_fit_model
    result <- biological_link_function(object, model_object, exact = FALSE)
    # Check that the result is a list with expected elements
    expect_true(is.list(result))
    expect_true("orthologs" %in% names(result))
    expect_true("singIST_samples" %in% names(result))
    expect_true("FC" %in% names(result))
})
