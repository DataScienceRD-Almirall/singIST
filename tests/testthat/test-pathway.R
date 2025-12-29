test_that("Check that slots of pathway class are consistent", {
    testthat::skip_on_cran()
    testthat::skip_on_bioc()
    # Error due to collection not being c2
    expect_error(create_pathway(standard_name = "KEGG_DENDRITIC",
                                dbsource = "KEGG",
                                collection = "x",
                                subcollection = "CP"))
    # Good case
    expect_list(create_pathway(standard_name = "KEGG_DENDRITIC",
                                dbsource = "KEGG",
                                collection = "c2",
                                subcollection = "CP"))
})
