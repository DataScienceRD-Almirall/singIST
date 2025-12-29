test_that("Check consistency of mapping.organism class slots", {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
    organism <- "Mus musculus"
    target_class <- "g1"
    base_class <- "g2"
    counts <- SeuratObject::pbmc_small # Toy dataset
    # Rename "orig.ident" variable to "donor"
    colnames(slot(counts, "meta.data"))[1] <- "donor"
    # Rename "group" variable to "class"
    colnames(slot(counts, "meta.data"))[6] <- "class"
    # Example existing mapping for T-cell but no mapping for Dendritic Cell
    celltype_mapping <- list("T-cell" = c("DETC", "T"),
    "Dendritic Cell" = c())
    # Rename "RNA_snn.res.1" variable to "celltype_cluster"
    colnames(slot(counts, "meta.data"))[7] <- "celltype_cluster"

    # Expect pass
    expect_list(create_mapping_organism(organism, target_class,
                                        base_class, celltype_mapping, counts))
    # Expect error due to inexisting "class", "celltype_cluster" and "donor"
    # variables
    colnames(slot(counts, "meta.data"))[6] <- "group"
    colnames(slot(counts, "meta.data"))[7] <- "celltype"
    colnames(slot(counts, "meta.data"))[1] <- "orig.ident"
    expect_error(create_mapping_organism(organism, target_class,
                                         base_class, celltype_mapping, counts))
})
