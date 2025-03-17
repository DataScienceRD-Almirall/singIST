test_that("Check that slots of hyperparameters class are consistent", {
    quantile_comb_table <- base::as.matrix(
        RcppAlgos::permuteGeneral(seq(0.05, 0.95, by = 0.50)),
        ncol = length(c("T-cell", "Dendritic Cell"))
    )
    # Error due to invalid outcome_type
    expect_error(new("hyperparameters",
                   quantile_comb_table = quantile_comb_table,
                   outcome_type = "notthisone",
                   number_PLS = as.integer(2),
                   folds_CV = as.integer(3),
                   repetition_CV = as.integer(1)))
    # Pass
    expect_class(new("hyperparameters",
                     quantile_comb_table = quantile_comb_table,
                     outcome_type = "multiclass",
                     number_PLS = as.integer(2),
                     folds_CV = as.integer(3),
                     repetition_CV = as.integer(1)),
                 "hyperparameters")
})
