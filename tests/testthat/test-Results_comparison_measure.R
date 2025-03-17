test_that("Results_comparison_measure calculates performance metrics correctly", {

    # Binary classification case
    Y_predict_binary <- c(1, 0, 1, 0, 1)
    Y_true_binary <- c(0, 0, 1, 1, 1)

    result_binary <- Results_comparison_measure(Y_predict_binary, Y_true_binary,
                                                outcome.type = "binary")

    # Check that the result is a named vector with the correct names
    expect_named(result_binary, c("accuracy", "balanced_accuracy", "precision",
                                    "recall", "F1"))
    # Test for multiclass classification
    Y_predict_multiclass <- c(1, 2, 2, 1, 3)
    Y_true_multiclass <- c(1, 2, 2, 3, 3)
    result_multiclass <- Results_comparison_measure(
        Y_predict_multiclass, Y_true_multiclass, outcome.type = "multiclass")
    # Check that the result is a named vector with the correct names
    expect_named(result_multiclass, c("accuracy", "balanced_accuracy",
                                        "precision", "recall", "F1"))
    # Test for edge case where all predictions are true positives
    # (perfect accuracy)
    Y_predict_perfect <- c(1, 1, 1, 1, 1)
    Y_true_perfect <- c(1, 1, 1, 1, 1)
    result_perfect <- Results_comparison_measure(
        Y_predict_perfect, Y_true_perfect, outcome.type = "binary")
    expect_equal(result_perfect["accuracy"][[1]], 1)
    expect_equal(result_perfect["precision"][[1]], 1)
    expect_equal(result_perfect["recall"][[1]], 1)
    expect_equal(result_perfect["F1"][[1]], 1)
    # Test for edge case where all predictions are false positives
    Y_predict_false <- c(0, 0, 0, 0, 0)
    Y_true_false <- c(1, 1, 1, 1, 1)
    result_false <- Results_comparison_measure(
        Y_predict_false, Y_true_false, outcome.type = "binary")
    expect_equal(result_false["accuracy"][[1]], 0)
    expect_equal(result_false["precision"][[1]], NaN)
    expect_equal(result_false["recall"][[1]], 0)
    expect_equal(result_false["F1"][[1]], NaN)
})
