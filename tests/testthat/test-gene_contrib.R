test_that("gene_contrib works as expected", {
    testthat::skip_on_cran()
    testthat::skip_if_not(interactive())
    file <- system.file("extdata", "example_superpathway_fit_model.rda", package = "singIST")
    load(file)
    file <- system.file("extdata", "example_mapping_organism.rda", package = "singIST")
    load(file)
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
    result <- gene_contrib(model, original$gene_contribution,
                            derived$gene_contribution, celltype)
    # Checking if the result is a data frame
    expect_true(is.data.frame(result))
    # Checking if the required columns are present
    expect_true("pathway" %in% colnames(result))
    expect_true("celltype" %in% colnames(result))
    expect_true("gene" %in% colnames(result))
    expect_true("contribution" %in% colnames(result))
})
