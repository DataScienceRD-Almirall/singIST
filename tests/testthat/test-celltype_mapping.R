test_that("celltype_mapping correctly updates celltype_cluster", {
    data("example_mapping_organism", package = "singIST")
    object <- example_mapping_organism
    # Ensure original celltype_cluster is updated
    new_object <- celltype_mapping(object)
    expect_true(all(names(new_object@celltype_mapping) %in%
                    new_object@counts$celltype_cluster))
    expect_false(any(is.na(new_object@counts$celltype_cluster)))
})
