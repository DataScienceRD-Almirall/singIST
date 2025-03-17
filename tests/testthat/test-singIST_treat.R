test_that("singIST_treat computes treated samples correctly", {
    data("example_mapping_organism", package = "singIST")
    object <- example_mapping_organism
    model_object <- example_superpathway_fit_model
    orthologs <- orthology_mapping(object, model_object)
    data <- celltype_mapping(object)
    slot(data, "counts")$test <- paste0(slot(data, "counts")$celltype_cluster,
                                        "_", slot(data, "counts")$class)
    SeuratObject::Idents(slot(data, "counts")) <- "test"
    logFC <- diff_expressed(data, exact = FALSE)
    singIST_samples <- singIST_treat(object, model_object, orthologs, logFC)
    # Check if output contains the correct components
    expect_true("singIST_samples" %in% names(singIST_samples))
    expect_true("FC" %in% names(singIST_samples))
    expect_true(is.list(singIST_samples$FC))
})
