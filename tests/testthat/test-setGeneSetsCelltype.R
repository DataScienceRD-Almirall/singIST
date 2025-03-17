test_that("Check that setGeneSetsCelltype works as intended", {
    data("example_superpathway_input", package = "singIST")
    # Throws warning of "Assuming organism to be human"
    expect_warning(result <-
                setRepeatGeneSets(example_superpathway_input@superpathway_info))
    ncelltypes <- length(example_superpathway_input@superpathway_info@celltypes)
    expect_true(all(lengths(result@gene_sets_celltype) ==
                           lengths(result@gene_sets_celltype)[1]))
    expect_true(length(result@gene_sets_celltype) == ncelltypes)
})
