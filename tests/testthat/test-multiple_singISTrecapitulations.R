test_that("Test multiple_singISTrecapitulations function", {
    data("example_superpathway_fit_model", package = "singIST")
    data("example_mapping_organism", package = "singIST")
    # Define example superpathway.input objects
    object1 <- object2 <- example_superpathway_input
    models <- list(object1, object2)
    # Test when using default arguments
    result_models <- multiple_fitOptimal(models)
    # Test when using default arguments for from_species
    result <- multiple_singISTrecapitulations(object, result_models,
                                                exact = FALSE)
    expect_equal(names(result), c("superpathway", "celltype", "gene", "FC"))
    expect_true(inherits(result$superpathway, "data.frame"))
    expect_true(inherits(result$celltype, "data.frame"))
    expect_true(inherits(result$gene, "data.frame"))
    # Test with custom from_species
    result <- multiple_singISTrecapitulations(
        object, result_models, from_species = list("hsapiens", "mmusculus"),
        exact = FALSE)
    expect_equal(length(result$superpathway$recapitulation), 2)
    # Test when the number of from_species does not match the number of models
    expect_error(multiple_singISTrecapitulations(
        object, model_objects, exact = FALSE, from_species = list("hsapiens")))
})
