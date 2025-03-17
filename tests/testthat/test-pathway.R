test_that("Check that slots of pathway class are consistent", {
    # Error due to collection not being c2
    expect_error(new("pathway",
                    standard_name = "KEGG_DENDRITIC",
                    dbsource = "KEGG",
                    collection = "x",
                    subcollection = "CP"
    ))
    # Error due to standard_name not containing dbsource
    expect_error(new("pathway",
                     standard_name = "DENDRITIC",
                     dbsource = "KEGG",
                     collection = "c2",
                     subcollection = "CP"
    ))
    # Good case
    expect_class(new("pathway",
                     standard_name = "KEGG_DENDRITIC",
                     dbsource = "KEGG",
                     collection = "c2",
                     subcollection = "CP"
    ), "pathway")
})
