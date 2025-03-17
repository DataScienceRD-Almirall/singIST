test_that("singISTrecapitulations works as expected", {
    data("example_superpathway_fit_model", package = "singIST")
    data("example_mapping_organism", package = "singIST")
    object <- example_mapping_organism
    model_object <- example_superpathway_fit_model
    result <- singISTrecapitulations(object, model_object, exact = FALSE)
    # Checking if the result is a list with correct components
    expect_true(is.list(result))
    expect_true("superpathway" %in% names(result))
    expect_true("celltype" %in% names(result))
    expect_true("gene" %in% names(result))
    expect_true("FC" %in% names(result))
    # Checking that the data frames have the expected columns
    expect_true("pathway" %in% colnames(result$superpathway))
    expect_true("recapitulation" %in% colnames(result$superpathway))
    expect_true("p_val" %in% colnames(result$superpathway))
    expect_true("celltype" %in% colnames(result$celltype))
    expect_true("recapitulation" %in% colnames(result$celltype))
    expect_true("orthology" %in% colnames(result$celltype))
    expect_true("pathway" %in% colnames(result$gene))
    expect_true("celltype" %in% colnames(result$gene))
    expect_true("gene" %in% colnames(result$gene))
    expect_true("contribution" %in% colnames(result$gene))
})
