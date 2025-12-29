test_that("singIST_treat computes treated samples correctly", {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
    file <- system.file("extdata", "example_mapping_organism.rda", package = "singIST")
    load(file)
    file <- system.file("extdata", "example_superpathway_fit_model.rda", package = "singIST")
    load(file)
    object <- example_mapping_organism
    model_object <- example_superpathway_fit_model
    orthologs <- orthology_mapping(object, model_object, from_species = "hsapiens")
    data <- celltype_mapping(object)
    data$counts$test <- paste0(data$counts$celltype_cluster,
                                        "_", data$counts$class)
    SeuratObject::Idents(data$counts) <- "test"
    logFC <- diff_expressed(data, exact = FALSE)
    singIST_samples <- singIST_treat(object, model_object, orthologs, logFC)
    # Check if output contains the correct components
    expect_true("singIST_samples" %in% names(singIST_samples))
    expect_true("FC" %in% names(singIST_samples))
    expect_true(is.list(singIST_samples$FC))
})
