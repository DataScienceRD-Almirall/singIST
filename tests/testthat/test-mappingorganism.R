test_that("Check consistency of mapping.organism class slots", {
    organism <- "Mus musculus"
    target_class <- "g1"
    base_class <- "g2"
    counts <- SeuratObject::pbmc_small # Toy dataset
    # Rename "group" variable to "class"
    colnames(slot(counts, "meta.data"))[6] <- "class"
    # Example existing mapping for T-cell but no mapping for Dendritic Cell
    celltype_mapping <- list("T-cell" = c("dγdT", "T"),
    "Dendritic Cell" = c())
    # Rename "RNA_snn.res.1" variable to "celltype_cluster"
    colnames(slot(counts, "meta.data"))[7] <- "celltype_cluster"

    # Expect pass
    expect_class(new("mapping.organism",
                     organism = organism,
                     target_class = target_class,
                     base_class = base_class,
                     celltype_mapping = celltype_mapping,
                     counts = counts), "mapping.organism")
    # Expect error due to inexisting "class" and "celltype_cluster" variables
    colnames(slot(counts, "meta.data"))[6] <- "group"
    colnames(slot(counts, "meta.data"))[7] <- "celltype"
    expect_error(new("mapping.organism",
                        organism = organism,
                        target_class = target_class,
                        base_class = base_class,
                        celltype_mapping = celltype_mapping,
                        counts = counts))
})
