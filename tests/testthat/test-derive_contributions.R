test_that("derive_contributions computes contributions correctly", {
    data("example_superpathway_fit_model", package = "singIST")
    model_object <- example_superpathway_fit_model
    data <- biological_link_function(
        example_mapping_organism,example_superpathway_fit_model,
        exact = FALSE)$singIST_samples
    # Call the function
    contributions_output <- derive_contributions(model_object, data)
    # Check if the output is a list
    expect_true(is.list(contributions_output))
    # Check the structure of the output
    expect_named(contributions_output, c("superpathway_score",
                                "celltype_contribution", "gene_contribution"))
    # Check if superpathway_score is a matrix
    expect_class(contributions_output$superpathway_score, "matrix")
    expect_equal(ncol(contributions_output$superpathway_score), nrow(data))
    # Check if celltype_contribution is a matrix
    expect_class(contributions_output$celltype_contribution, "matrix")
    expect_equal(nrow(contributions_output$celltype_contribution),
            length(model_object@superpathway_input@superpathway_info@celltypes))
    # Check if gene_contribution is a list of matrices
    expect_type(contributions_output$gene_contribution, "list")
    expect_true(all(sapply(contributions_output$gene_contribution, is.matrix)))
})

test_that("derive_contributions handles edge cases", {
    data("example_superpathway_fit_model", package = "singIST")
    # Test with an empty data frame or NA values in the data
    empty_data <- data.frame(matrix(ncol = 0, nrow = 0))
    expect_error(derive_contributions(example_superpathway_fit_model,
                                        empty_data))
})

