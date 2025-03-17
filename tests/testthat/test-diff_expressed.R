test_that("diff_expressed computes differential expression", {
    data("example_mapping_organism", package = "singIST")
    object <- example_mapping_organism
    data <- celltype_mapping(object)
    slot(data, "counts")$test <- paste0(slot(data, "counts")$celltype_cluster,
    "_", slot(data, "counts")$class)
    SeuratObject::Idents(slot(data, "counts")) <- "test"
    logFC <- diff_expressed(data, exact = FALSE)
    # Check if output is a list and contains the expected elements
    expect_true(is.list(logFC))
    expect_true(length(logFC) > 0)
    expect_true(all(sapply(logFC, is.data.frame)))
})
