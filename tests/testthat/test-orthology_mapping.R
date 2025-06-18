test_that("orthology_mapping produces correct output", {
    file <- system.file("extdata", "example_mapping_organism.rda", package = "singIST")
    load(file)
    file <- system.file("extdata", "example_superpathway_fit_model.rda", package = "singIST")
    load(file)
    object <- example_mapping_organism
    model_object <- example_superpathway_fit_model
    orthologs <- orthology_mapping(object, model_object, from_species = "hsapiens")
    # Check if output is a list
    expect_true(is.list(orthologs))
    expect_true(length(orthologs) > 0)
})

