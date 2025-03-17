test_that("gene_contrib works as expected", {
    data("example_superpathway_fit_model", package = "singIST")
    data("example_mapping_organism", package = "singIST")
    model <- example_superpathway_fit_model
    mapped <- example_mapping_organism
    singIST_samples <- biological_link_function(mapped, model,
                                                exact = FALSE)$singIST_samples
    original <- derive_contributions(model, singIST_samples)
    derived <- derive_contributions(model,
                                    slot(model, "model_fit")$predictor_block)
    celltype <- celltype_recap(model, original$celltype_contribution,
                             derived$celltype_contribution)
    # Running the function
    result <- gene_contrib(model_object, original$gene_contribution,
                            derived$gene_contribution, celltype)
    # Checking if the result is a data frame
    expect_true(is.data.frame(result))
    # Checking if the required columns are present
    expect_true("pathway" %in% colnames(result))
    expect_true("celltype" %in% colnames(result))
    expect_true("gene" %in% colnames(result))
    expect_true("contribution" %in% colnames(result))
})
