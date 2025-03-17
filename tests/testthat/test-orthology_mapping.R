test_that("orthology_mapping produces correct output", {
    data("example_mapping_organism", package = "singIST")
    object <- example_mapping_organism
    model_object <- example_superpathway_fit_model
    orthologs <- orthology_mapping(object, model_object)
    # Check if output is a list
    expect_true(is.list(orthologs))
    expect_true(length(orthologs) > 0)
})

