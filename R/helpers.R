#' @title Update block of predictor matrices in matrixToBlock()
#'
#' @description
#' Fill up matrix with the corresponding expression values
#'
#' @param celltype Cell types modelled
#' @param observed_gene_sets Gene sets observed from your dataset
#' @param block_predictor Block of predictor matrices to update
#' @param matrix To iteratively update with block_predictor values
#' @name helpers
#' @rdname helpers
#' @export
update_block <- function(celltype, observed_gene_sets,
                            block_predictor = block_predictor,
                            matrix = matrix){
    matching_rows <- base::grep(paste0("^", celltype, "_"),
                                rownames(matrix), value = TRUE)
    filtered_matrix <- matrix[matching_rows,
                                observed_gene_sets, drop = FALSE]
    rownames(filtered_matrix) <- base::gsub(paste0("^", celltype, "_"),
                                            "", matching_rows)
    colnames(filtered_matrix) <- paste0(celltype, "*", observed_gene_sets)
    common_columns <- intersect(colnames(filtered_matrix),
                                colnames(block_predictor))
    block_predictor[, common_columns] <- filtered_matrix[, common_columns]
    return(block_predictor)
}

#' @title Get index measure
#' @description
#' \code{get_measure_index()} returns the index associated to each performance
#' measure
#' @param measure A character either "accuracy", "B_accuracy", "precision",
#' "recall" or "F1.
#' @name helpers
#' @rdname helpers
#' @export
#' @examples
#' measure <- "F1"
#' get_measure_index(measure)
get_measure_index <- function(measure){
    if(measure == "accuracy") {
        measure_selected <- 1
    }
    if(measure == "B_accuracy") {
        measure_selected <- 2
    }
    if(measure == "precision") {
        measure_selected <- 3
    }
    if(measure == "recall") {
        measure_selected <- 4
    }
    if(measure == "F1") {
        measure_selected <- 5
    }
    return(measure_selected)
}

#' @title Get Training and Validation Sets
#'
#' @description Splits the predictor and response matrices into training and
#' validation sets for leave-one-out cross-validation.
#'
#' @param X.matrix Predictor matrix.
#' @param Y.matrix Response matrix.
#' @param validation_index Index of the validation sample.
#'
#' @return A list containing the training and validation sets:
#' \item{E_matrix_validation}{Validation predictor matrix}
#' \item{F_matrix_validation}{Validation response matrix}
#' \item{E_matrix_training}{Training predictor matrix}
#' \item{F_matrix_training}{Training response matrix}
#' @name helpers
#' @rdname helpers
#' @export
#' @examples
#' X <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' Y <- matrix(sample(0:1, 10, replace = TRUE), ncol = 1)
#' result <- get_train_val_sets(X, Y, validation_index = 2)
#' str(result)
get_train_val_sets <- function(X.matrix, Y.matrix, validation_index) {
    aux <- seq_len(nrow(Y.matrix))
    training_index <- aux[aux != validation_index]
    return(list(
        "E_matrix_validation" = X.matrix[validation_index, , drop = FALSE],
        "F_matrix_validation" = Y.matrix[validation_index, , drop = FALSE],
        "E_matrix_training" = X.matrix[training_index, , drop = FALSE],
        "F_matrix_training" = Y.matrix[training_index, , drop = FALSE]
    ))
}

#' @title Parallelize quantile hyperparameter tuning in LOOCV
#' @description
#' Function to train and validate asmbPLSDA excluding one observation
#' parallelized for each quantile combination provided
#' @param j Observation to exclude from training
#' @param ... Other parameters to pass onto `furrr::future_pmap`
#' @param results_CV_summary_n Matrix to store LOOCV results
#' @param F_matrix_validation_bind Matrix to validate LOOCV iteration
#' @param X.matrix A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param Y.matrix A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param PLS_term A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param X.dim A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param quantile.comb.table A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param outcome.type A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param quantile_table_CV A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param K A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param n_quantile_comb A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param Method A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param measure A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param expected.measure.increase A parameter passed from
#' \link{asmbPLSDA.cv.loo}
#' @param center A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param scale A parameter passed from \link{asmbPLSDA.cv.loo}
#' @param maxiter A parameter passed from \link{asmbPLSDA.cv.loo}
#' @name helpers
#' @rdname helpers
#'
#' @import asmbPLS
#' @export
quantile_computation <-
    function(j, ..., results_CV_summary_n, F_matrix_validation_bind,
                X.matrix, Y.matrix, PLS_term = 1, X.dim, quantile.comb.table,
                outcome.type = c("binary", "multiclass"), quantile_table_CV, K,
                n_quantile_comb, Method = NULL, measure = "B_accuracy",
                expected.measure.increase = 0.005, center = TRUE, scale = TRUE,
                maxiter = 100){
            validation_index <- j
            aux <- seq_len(K)
            training_index <- aux[aux != j]
            E_matrix_validation <- X.matrix[validation_index, , drop = FALSE]
            F_matrix_validation <- Y.matrix[validation_index, , drop = FALSE]
            E_matrix_training <- X.matrix[training_index, , drop = FALSE]
            F_matrix_training <- Y.matrix[training_index, , drop = FALSE]
            # calculate overall/balanced accuracy using different quantile
            # combinations
            for (l in seq_len(n_quantile_comb)) {
                quantile_table_CV[PLS_term, seq_len(length(X.dim))] <-
                    quantile.comb.table[l, seq_len(length(X.dim))]
                quantile_temp <-
                    quantile_table_CV[seq_len(PLS_term), seq_len(length(X.dim)),
                                        drop = FALSE]
                asmbPLSDA_fit_results <-
                    asmbPLS::asmbPLSDA.fit(
                        E_matrix_training,
                        F_matrix_training, PLS_term, X.dim,
                        quantile_temp, outcome.type, center, scale, maxiter)
                asmbPLSDA_predict_results <-
                    asmbPLS::asmbPLSDA.predict(
                        asmbPLSDA_fit_results, E_matrix_validation,
                        PLS_term, Method)
                Y_pred <- as.numeric(asmbPLSDA_predict_results["Y_pred"])
                results_CV_summary_n[l, j] <- Y_pred
                F_matrix_validation_bind[l, j] <- F_matrix_validation
            }
            return(list("results_CV_summary_n" = results_CV_summary_n,
                        "F_matrix_validation_bind" = F_matrix_validation_bind,
                        "obs" = j)
            )
}

#' @title Evaluate Quantile Combinations
#'
#' @description Computes the prediction accuracy for different quantile
#' combinations by fitting the asmbPLSDA model and making predictions.
#'
#' @param j PLS iteration from \code{execute_sequential_cv()}
#' @param results_CV_summary_n Matrix to store the predicted LOOCV samples for
#' each quantile combination
#' @param F_matrix_validation_bind Matrix to store the validation LOOCV samples
#' for each quantile combination
#' @param E_matrix_training Training predictor matrix.
#' @param F_matrix_training Training response matrix.
#' @param E_matrix_validation Validation predictor matrix.
#' @param F_matrix_validation Validation response matrix
#' @param quantile_table_CV Table storing quantile combinations.
#' @param i Number of PLS components.
#' @param X.dim List of gene set sizes for each cell type.
#' @param quantile.comb.table Table of quantile combinations.
#' @param outcome.type Character indicating outcome type (`"binary"` or
#' `"multiclass"`).
#' @param center Logical; whether to center data.
#' @param scale Logical; whether to scale data.
#' @param maxiter Integer; maximum number of iterations.
#' @param Method Prediction method.
#' @name helpers
#' @rdname helpers
#' @return A numeric vector containing predicted values for validation samples.
#' @import asmbPLS
#' @export
#'
#' @examples
#' E_train <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' F_train <- matrix(sample(0:1, 10, replace = TRUE), ncol = 1)
#' E_valid <- matrix(rnorm(10), nrow = 1, ncol = 10)
#' F_valid <- matrix(1, nrow = 1, ncol = 1)
#' quantile_table <- matrix(runif(2), nrow = 1, ncol = 2)
#' quantile_table_CV <- matrix(runif(7), nrow = 1, ncol = 7)
#' results_CV_summary_n <- matrix(0, nrow = 1, ncol = 2)
#' F_matrix_validation_bind <- matrix(0, nrow = 1, ncol = 2)
#' result <- evaluate_quantile_combinations(j=1, E_matrix_training = E_train,
#'                                          F_matrix_training = F_train,
#'                                          E_matrix_validation = E_valid,
#'                                          F_matrix_validation = F_valid,
#'                                          F_matrix_validation_bind =
#'                                          F_matrix_validation_bind,
#'                                          results_CV_summary_n =
#'                                          results_CV_summary_n,
#'                                          quantile_table_CV=quantile_table_CV,
#'                                          i = 1, X.dim = c(5,5),
#'                                          quantile.comb.table =quantile_table,
#'                                          outcome.type = "binary",
#'                                          center = TRUE,
#'                                          scale = TRUE, maxiter = 100,
#'                                          Method = NULL)
#' print(result)
evaluate_quantile_combinations <- function(j, results_CV_summary_n,
                                            F_matrix_validation_bind,
                                            E_matrix_training,
                                            F_matrix_training,
                                            E_matrix_validation,
                                            F_matrix_validation,
                                            quantile_table_CV, i, X.dim,
                                            quantile.comb.table, outcome.type,
                                            center, scale, maxiter, Method) {
        for(l in seq_len(nrow(quantile.comb.table))){
            quantile_table_CV[i, seq_len(length(X.dim))] <-
                quantile.comb.table[l, seq_len(length(X.dim))]
            quantile_temp <- quantile_table_CV[
                seq_len(i), seq_len(length(X.dim)), drop = FALSE]
            fit_results <- asmbPLS::asmbPLSDA.fit(
                E_matrix_training, F_matrix_training, i, X.dim, quantile_temp,
                outcome.type, center, scale,maxiter)
            predict_results <- asmbPLS::asmbPLSDA.predict(fit_results,
                                                            E_matrix_validation,
                                                            i, Method)
            Y_pred <- as.numeric(predict_results["Y_pred"])
            results_CV_summary_n[l, j] <- Y_pred
            F_matrix_validation_bind[l, j] <- F_matrix_validation
        }
    return(list("results_CV_summary_n" = results_CV_summary_n,
                "F_matrix_validation_bind" = F_matrix_validation_bind))
}

#' @title Execute Parallel Cross-Validation
#'
#' @description Performs leave-one-out cross-validation (LOO-CV) in parallel.
#'
#' @param K Number of samples.
#' @param results_CV_summary_n Matrix to store CV results.
#' @param F_matrix_validation_bind Matrix to store validation responses.
#' @param X.matrix Predictor matrix.
#' @param Y.matrix Response matrix.
#' @param PLS_term Number of PLS components.
#' @param X.dim List of gene set sizes.
#' @param quantile.comb.table Table of quantile combinations.
#' @param outcome.type Outcome type (`"binary"` or `"multiclass"`).
#' @param quantile_table_CV Table storing quantile results.
#' @param Method Prediction method.
#' @param measure Performance metric.
#' @param expected.measure.increase Minimum accuracy increase threshold.
#' @param center Logical; whether to center data.
#' @param scale Logical; whether to scale data.
#' @param maxiter Maximum number of iterations.
#' @param Method Method to predict class
#' @param BPPARAM A `BiocParallel::bpparam()` with parallelization options
#' @name helpers
#' @rdname helpers
#' @return A list containing updated `results_CV_summary_n` and
#' `F_matrix_validation_bind` matrices.
#' @import BiocParallel
#' @export
#' @examples
#' set.seed(123)
#' K <- 5
#' X <- matrix(rnorm(50), nrow = 5, ncol = 10)
#' Y <- matrix(sample(0:1, 5, replace = TRUE), ncol = 1)
#' quantile_comb_table <- matrix(runif(10), nrow = 2, ncol = 10)
#' results_CV_summary_n <- matrix(0, nrow = 2, ncol = K)
#' F_matrix_validation_bind <- matrix(0, nrow = 2, ncol = K)
#' # Parallelization options
#' library(BiocParallel)
#' register(SnowParam(workers = 2, exportglobals = FALSE, progressbar = TRUE),
#' default = TRUE)
#' output <- execute_parallel_cv(K, results_CV_summary_n,
#'                               F_matrix_validation_bind, X, Y, PLS_term = 1,
#'                               X.dim = c(5,5),
#'                               quantile.comb.table = quantile_comb_table,
#'                               outcome.type = "binary",
#'                               quantile_table_CV = quantile_comb_table,
#'                               measure = "B_accuracy",
#'                               expected.measure.increase = 0.005,
#'                               center = TRUE, scale = TRUE, maxiter = 100,
#'                               Method = NULL)
#'register(SerialParam(), default = TRUE) # disable parallelization
#' str(output)
execute_parallel_cv <- function(K, results_CV_summary_n,
                                F_matrix_validation_bind, X.matrix, Y.matrix,
                                PLS_term, X.dim, quantile.comb.table,
                                outcome.type, quantile_table_CV, Method,
                                measure, expected.measure.increase, center,
                                scale, maxiter,
                                BPPARAM = BiocParallel::bpparam()) {
    output <- BiocParallel::bplapply(
        seq_len(K), quantile_computation,
        results_CV_summary_n = results_CV_summary_n,
        F_matrix_validation_bind = F_matrix_validation_bind,X.matrix = X.matrix,
        Y.matrix = Y.matrix, PLS_term = PLS_term, X.dim = X.dim,
        quantile.comb.table = quantile.comb.table, outcome.type = outcome.type,
        quantile_table_CV = quantile_table_CV, K = K,
        n_quantile_comb = nrow(quantile.comb.table), Method = Method,
        measure = measure,expected.measure.increase = expected.measure.increase,
        center = center, scale = scale, maxiter = maxiter,
        BPPARAM = BPPARAM
    )
    results_CV_summary_n <- output[[1]]$results_CV_summary_n
    F_matrix_validation_bind <- output[[1]]$F_matrix_validation_bind
    for (ncol in 2:K) {
        results_CV_summary_n[, ncol]<-output[[ncol]]$results_CV_summary_n[,ncol]
        F_matrix_validation_bind[, ncol] <-
            output[[ncol]]$F_matrix_validation_bind[, ncol]
    }
    return(list("results_CV_summary_n" = results_CV_summary_n,
                "F_matrix_validation_bind" = F_matrix_validation_bind))
}

#' @title Compute performance measures for one PLS and
#' each quantile combination
#' @description
#' Computes the performance measure selected between the training LOOCV samples
#' and the validation LOOCV samples for all the quantile combination
#' @param n_quantile_comb Passed from \link{asmbPLSDA.cv.loo}
#' @param results_CV_summary_n Passed from \link{asmbPLSDA.cv.loo}
#' @param F_matrix_validation_bind Passed from \link{asmbPLSDA.cv.loo}
#' @param outcome.type Passed from \link{asmbPLSDA.cv.loo}
#' @param measure_selected Passed from \link{asmbPLSDA.cv.loo}
#' @name helpers
#' @rdname helpers
#' @returns
#' A vector with the performance measure of each quantile combination
#' @export
performance_measures <- function(n_quantile_comb, results_CV_summary_n,
                                    F_matrix_validation_bind, outcome.type,
                                    measure_selected){
    measure_acc <- c()
    for(l in seq_len(n_quantile_comb)){
        Y_pred <- as.vector(results_CV_summary_n[l, ])
        F_matrix_validation <- as.vector(F_matrix_validation_bind[l, ])
        measure_new <- Results_comparison_measure(Y_pred, F_matrix_validation,
                                                    outcome.type)
        measure_acc <- c(measure_acc, measure_new[measure_selected])
    }
    return(measure_acc)
}

#' @title Compute performance metrics with the optimal quantiles and PLS
#' components
#' @description
#' For an optimal quantile combination and PLS component it computes its
#' performance metrics between the training and validation sets
#'
#' @param K Passed from \link{asmbPLSDA.cv.loo}
#' @param X.matrix Passed from \link{asmbPLSDA.cv.loo}
#' @param Y.matrix Passed from \link{asmbPLSDA.cv.loo}
#' @param i Passed from \link{asmbPLSDA.cv.loo}
#' @param X.dim Passed from \link{asmbPLSDA.cv.loo}
#' @param quantile_table_CV Passed from \link{asmbPLSDA.cv.loo}
#' @param outcome.type Passed from \link{asmbPLSDA.cv.loo}
#' @param center Passed from \link{asmbPLSDA.cv.loo}
#' @param scale Passed from \link{asmbPLSDA.cv.loo}
#' @param maxiter Passed from \link{asmbPLSDA.cv.loo}
#' @param Method Passed from \link{asmbPLSDA.cv.loo}
#' @name helpers
#' @rdname helpers
#' @import asmbPLS
#' @returns
#' Optimal quantile table for each PLS with all its performance measures
#' @export
compute_final_measures <- function(K, X.matrix, Y.matrix, i, X.dim,
                                    quantile_table_CV, outcome.type,
                                    center, scale, maxiter, Method){
    Y_pred_bind <- matrix()
    F_matrix_validation_bind <- matrix()
    for(j in seq_len(K)){
        train_val <- get_train_val_sets(X.matrix, Y.matrix, j)
        quantile_temp <- quantile_table_CV[seq_len(i), seq_len(length(X.dim)),
                                            drop = FALSE]
        # Fit model using the training set
        asmbPLSDA_fit_results <- asmbPLS::asmbPLSDA.fit(
            train_val$E_matrix_training, train_val$F_matrix_training, i, X.dim,
            quantile_temp, outcome.type, center, scale, maxiter)
        asmbPLSDA_predict_results <- asmbPLS::asmbPLSDA.predict(
            asmbPLSDA_fit_results, train_val$E_matrix_validation, i, Method)
        Y_pred <- as.matrix(asmbPLSDA_predict_results["Y_pred"])
        rownames(Y_pred) <- NULL
        Y_pred_bind <- rbind(Y_pred_bind, Y_pred)
        F_matrix_validation_bind <- rbind(F_matrix_validation_bind,
                                            train_val$F_matrix_validation)
    }
    # Avoid the first row which is always NA
    Y_pred_bind <- Y_pred_bind[2:nrow(Y_pred_bind), , drop = FALSE]
    F_matrix_validation_bind <- F_matrix_validation_bind[
        2:nrow(F_matrix_validation_bind), , drop = FALSE]
    # Compute the performance metrics for the validation and training sets
    measure <- Results_comparison_measure(Y_pred_bind, F_matrix_validation_bind,
                                            outcome.type)
    quantile_table_CV[i, (length(X.dim)+1):ncol(quantile_table_CV)] <- measure
    # Update colnames of optimal quantile table including blocks
    colnames(quantile_table_CV)[(length(X.dim)+1):ncol(quantile_table_CV)] <-
        names(measure)
    colnames(quantile_table_CV)[seq_len(length(X.dim))] <-
        paste(rep("block", length(X.dim)), seq_len(length(X.dim)))
    return(quantile_table_CV)
}


#' @title Compute optimal number of PLS
#' @description
#' Selects the optimal number of PLS according to the performance measure
#'
#' @param PLS_term Passed from \link{asmbPLSDA.cv.loo}
#' @param quantile_table_CV Passed from \link{asmbPLSDA.cv.loo}
#' @param X.dim Passed from \link{asmbPLSDA.cv.loo}
#' @param measure_selected Passed from \link{asmbPLSDA.cv.loo}
#' @param expected.measure.increase Passed from \link{asmbPLSDA.cv.loo}
#' @name helpers
#' @rdname helpers
#' @returns
#' An integer with the optimal number of PLS
#'
#' @export
select_optimal_PLS <- function(PLS_term, quantile_table_CV, X.dim,
                                measure_selected, expected.measure.increase){
    optimal_nPLS <- 1
    if(PLS_term > 1){
        for(i in seq_len(PLS_term-1)){
            col <- length(X.dim)+measure_selected
            current_measure <- quantile_table_CV[i,col]
            next_measure <- quantile_table_CV[i+1, col]
            if(next_measure > current_measure+expected.measure.increase){
                optimal_nPLS <- optimal_nPLS + 1
            }else{break}
        }
    }
    return(optimal_nPLS)
}

#' @title Execute sequential iterations for LOOCV
#'
#' @description
#' Iterates over all quantiles to generate the fitted asmbPLSDA for each and
#' its associated predicted values
#'
#' @param K Passed from \link{asmbPLSDA.cv.loo}
#' @param n_quantile_comb Passed from \link{asmbPLSDA.cv.loo}
#' @param X.matrix Passed from \link{asmbPLSDA.cv.loo}
#' @param results_CV_summary_n Passed from \link{asmbPLSDA.cv.loo}
#' @param F_matrix_validation_bind Passed from \link{asmbPLSDA.cv.loo}
#' @param Y.matrix Passed from \link{asmbPLSDA.cv.loo}
#' @param PLS_term Passed from \link{asmbPLSDA.cv.loo}
#' @param X.dim Passed from \link{asmbPLSDA.cv.loo}
#' @param quantile.comb.table Passed from \link{asmbPLSDA.cv.loo}
#' @param outcome.type Passed from \link{asmbPLSDA.cv.loo}
#' @param quantile_table_CV Passed from \link{asmbPLSDA.cv.loo}
#' @param measure Passed from \link{asmbPLSDA.cv.loo}
#' @param expected.measure.increase Passed from \link{asmbPLSDA.cv.loo}
#' @param center Passed from \link{asmbPLSDA.cv.loo}
#' @param scale Passed from \link{asmbPLSDA.cv.loo}
#' @param maxiter Passed from \link{asmbPLSDA.cv.loo}
#' @param Method Passed from \link{asmbPLSDA.cv.loo}
#' @name helpers
#' @rdname helpers
#' @returns
#' A list with the true class of each LOOCV sample and its predicted class
#' for each quantile combination
#' @export
execute_sequential_cv <- function(
        K, n_quantile_comb, results_CV_summary_n, F_matrix_validation_bind,
        X.matrix, Y.matrix, PLS_term, X.dim, quantile.comb.table, outcome.type,
        quantile_table_CV, measure, expected.measure.increase, center, scale,
        maxiter, Method) {
    for (j in seq_len(K)) {
        train_val <- get_train_val_sets(X.matrix, Y.matrix, j)
        evaluation_quantiles <- evaluate_quantile_combinations(
                    j, results_CV_summary_n, F_matrix_validation_bind,
                    train_val$E_matrix_training,
                    train_val$F_matrix_training,
                    train_val$E_matrix_validation,
                    train_val$F_matrix_validation,
                    quantile_table_CV, PLS_term, X.dim,
                    quantile.comb.table, outcome.type, center, scale, maxiter,
                    Method
                    )
        results_CV_summary_n[ , j] <-
            evaluation_quantiles$results_CV_summary_n[, j]
        F_matrix_validation_bind[ , j] <-
            evaluation_quantiles$F_matrix_validation_bind[, j]
    }
    return(list("results_CV_summary_n" = results_CV_summary_n,
                "F_matrix_validation_bind" = F_matrix_validation_bind))
}

#' @title Initialize result storage for permutations
#'
#' @description Creates a structured list to store permutation results.
#'
#' @param npermut Number of permutations.
#' @param q Number of classes.
#'
#' @return A list containing initialized data frames for permutation statistics.
#' @name helpers
#' @rdname helpers
#' @export
#' @examples
#' initialize_results(100, 3)
initialize_results <- function(npermut, q) {
    dimlabP <- c("NoPermut", paste0("permut", seq_len(npermut)))
    res <- list(
        RV.YYpermut.values = data.frame(dimlabP, rep(NA, npermut + 1)),
        cor.YYpermut.values = data.frame(
            dimlabP, matrix(NA, ncol = q, nrow = npermut + 1)),
        prctGlob.Ychange.values = data.frame(
            dimlabP, rep(NA, npermut + 1)),
        prct.Ychange.values = data.frame(
            dimlabP,
            matrix(NA,ncol = q + 1,
            nrow = npermut + 1))
    )
    return(res)
}

#' @title Permute the response matrix Y
#'
#' @description Performs random permutations of the response matrix.
#'
#' @param Y.matrix Original response matrix.
#' @param nr Number of samples.
#' @param nbObsPermut Number of samples to permute per iteration.
#' @param j Current permutation index.
#' @name helpers
#' @rdname helpers
#' @return A permuted response matrix.
#' @export
#' @examples
#' permute_Y_matrix(matrix(rnorm(100), 10, 10), nr = 10, nbObsPermut = 3, j = 2)
permute_Y_matrix <- function(Y.matrix, nr, nbObsPermut, j) {
    if (j == 1) return(Y.matrix)
    nObsP <- if (is.null(nbObsPermut)) sample(nr, 1) else nbObsPermut
    for (o in seq_len(nObsP)) {
        inds <- sample(nr, 2)
        Y.matrix[inds, ] <- Y.matrix[rev(inds), ]
    }
    return(Y.matrix)
}

#' @title Compute permutation statistics
#'
#' @description Calculates correlation, percentage change, and RV coefficient.
#'
#' @param res List of results to store statistics.
#' @param Y.matrix Original response matrix.
#' @param Ypermut Permuted response matrix.
#' @param j Current permutation index.
#' @param q Number of classes.
#' @param nr Number of samples.
#'
#' @return Updated result list with permutation statistics.
#' @name helpers
#' @rdname helpers
#' @import FactoMineR stats
#' @export
#'
#' @examples
#' res <- initialize_results(100, 3)
#' compute_permutation_stats(res, matrix(rnorm(100), 10, 10),
#' matrix(rnorm(100), 10, 10), j = 2, q = 3, nr = 10)
compute_permutation_stats <- function(res, Y.matrix, Ypermut, j, q, nr) {
    res$cor.YYpermut.values[j, -1] <- base::vapply(
        seq_len(q),
        function(Q) stats::cor(Y.matrix[, Q], Ypermut[, Q]),
        FUN.VALUE = numeric(1))
    res$prct.Ychange.values[j, -1] <- base::colMeans(Y.matrix != Ypermut)
    res$RV.YYpermut.values[j, 2] <- FactoMineR::coeffRV(Y.matrix, Ypermut)$rv
    res$prctGlob.Ychange.values[j, 2] <- base::mean((Y.matrix - Ypermut)^2)
    return(res)
}

#' @title Randomly select samples for cross-validation
#'
#' @description Selects sample indices for training and validation.
#'
#' @param object The fitted model object.
#' @param nr Number of samples.
#' @param Nc Number of samples to drop at each permutation.
#' @name helpers
#' @rdname helpers
#' @return A vector of selected sample indices.
#' @export
select_samples <- function(object, nr, Nc) {
    if (Nc == 1) return(sample(nr, Nc))
    table_counts <- table(object@superpathway_input@sample_class)/nr
    prob <- as.vector(table_counts[object@superpathway_input@sample_class])
    s <- sample(nr, Nc, prob = prob)
    return(s)
}

#' @title Fit a model on permuted data
#'
#' @description Fits the asmbPLS-DA model using permuted data.
#'
#' @param object The fitted model object.
#' @param X_train Training predictor matrix.
#' @param Y_train Training response matrix.
#' @param maxiter Maximum number of iterations.
#' @name helpers
#' @rdname helpers
#' @return The fitted asmbPLS-DA model.
#' @import asmbPLS
#' @export
fit_permuted_model <- function(object, X_train, Y_train, maxiter) {
    model <- asmbPLS::asmbPLSDA.fit(
        X_train, Y_train, object@hyperparameters_fit@number_PLS,
        lengths(object@model_fit$observed_gene_sets),
        object@hyperparameters_fit@quantile_comb_table, center = TRUE,
        scale = TRUE, object@hyperparameters_fit@outcome_type, maxiter
    )
    return(model)
}

#' @title Compute performance of permuted asmbPLSDA model against ground truth
#'
#' @param res List of results to store statistics
#' @param Modelpermut Permuted asmbPLSDA model
#' @param X_train Training predictor blocks
#' @param X_val Validation predictor blocks
#' @param Y.matrix True response matrix
#' @param s Validation samples
#' @param measure Selected measure to compute performance
#' @param j Current permutation index
#' @param nr Number of samples
#' @param object A \link{superpathway.fit.model-class} object
#' @param Method Method to compute predictions
#' @name helpers
#' @rdname helpers
#' @returns
#' Res list including the performance measure of the permuted model
#' @export
evaluate_performance <- function(res, Modelpermut, X_train, X_val,
                                Y.matrix, s, measure, j, nr, Method, object) {
    Yperm_pred <- numeric(nr)
    Yperm_pred[s] <- asmbPLS::asmbPLSDA.predict(
        Modelpermut, X_val, object@hyperparameters_fit@number_PLS,
        Method)$Y_pred
    Yperm_pred[-s] <- asmbPLS::asmbPLSDA.predict(
        Modelpermut, X_train, object@hyperparameters_fit@number_PLS,
        Method)$Y_pred
    measures <- Results_comparison_measure(
        Yperm_pred, Y.matrix, object@hyperparameters_fit@outcome_type)
    res$prct.Ychange.values[j, ncol(res$prct.Ychange.values)] <-
        measures[get_measure_index(measure)]
    return(res)
}

#' @title Compute the p-value for permutation test
#'
#' @description Computes the p-value for the observed CV error against the null
#' distribution of errors generated from permutation testing.
#'
#' @param null_errors A vector of errors from the null distribution
#' (permuted errors).
#' @param CV_error The observed cross-validation error.
#'
#' @return The computed p-value.
#' @name helpers
#' @rdname helpers
#' @import stats
#' @export
#' @examples
#' null_errors <- c(0.3, 0.4, 0.35, 0.33)
#' CV_error <- 0.32
#' compute_pvalue(null_errors, CV_error)
compute_pvalue <- function(null_errors, CV_error) {
    ecdf_errors <- stats::ecdf(null_errors)
    pvalue <- 1 - ecdf_errors(CV_error)
    return(pvalue)
}

#' @title Compute the 95% confidence interval for null distribution
#'
#' @description Calculates the 95% confidence interval for the null distribution
#' of permutation errors.
#'
#' @param m A vector of errors from the null distribution (permuted errors).
#'
#' @return A numeric vector containing the lower and upper bounds of the 95%
#' confidence interval.
#' @name helpers
#' @rdname helpers
#' @import stats
#' @export
#' @examples
#' null_errors <- c(0.3, 0.4, 0.35, 0.33)
#' compute_IC95(null_errors)
compute_IC95 <- function(m) {
    rt <- c(
        round(stats::quantile(m, 0.05), 5),
        round(stats::quantile(m, 0.95), 5)
    )
    return(rt)
}

#' @title Jackknife Resampling for CIP/GIP Calculation
#' @description Perform the jackknife resampling procedure for CIP/GIP
#' calculations.
#' @param object Model object.
#' @param X.matrix Predictor matrix.
#' @param Y.matrix Response matrix.
#' @param K Number of samples.
#' @param maxiter Maximum iterations.
#' @param X.dim Dimensions of the observed gene sets.
#' @name helpers
#' @rdname helpers
#' @export
#' @returns A list with the observed CIP and GIP distributions.
jackknife_CIP_GIP <- function(object, X.matrix, Y.matrix, K, maxiter, X.dim) {
    CIP_GIP_variability <- list()
    for (j in seq_len(K)) {
        # Resample the training set
        training_index <- setdiff(seq_len(K), j)
        E_matrix_training <- X.matrix[training_index, , drop = FALSE]
        F_matrix_training <- Y.matrix[training_index, , drop = FALSE]
        # Fit model
        asmbPLSDA_fit_results <- asmbPLS::asmbPLSDA.fit(
            X.matrix = E_matrix_training, Y.matrix = F_matrix_training,
            PLS.comp = object@hyperparameters_fit@number_PLS, X.dim = X.dim,
            quantile.comb = object@hyperparameters_fit@quantile_comb_table,
            outcome.type = object@hyperparameters_fit@outcome_type,
            center = TRUE, scale = TRUE, maxiter = maxiter)
        object@model_fit$`asmbPLS-DA` <- asmbPLSDA_fit_results
        aux <- CIP_GIP(object)
        # Store results
        if (j == 1) {
            CIP_GIP_variability <- aux
        } else {
            for (b in seq_len(length(X.dim))) {
                CIP_GIP_variability$GIP[[b]] <- cbind(
                    CIP_GIP_variability$GIP[[b]], aux$GIP[[b]])
            }
            CIP_GIP_variability$CIP <- cbind(CIP_GIP_variability$CIP, aux$CIP)
        }
    }
    return(CIP_GIP_variability)
}

#' @title Subsampling Resampling for CIP/GIP Calculation
#' @description Perform the subsampling procedure for CIP/GIP calculations.
#' @param object Model object.
#' @param X.matrix Predictor matrix.
#' @param Y.matrix Response matrix.
#' @param K Number of samples.
#' @param M Number of classes.
#' @param nsubsampling Number of subsamples.
#' @param maxiter Maximum iterations.
#' @param X.dim Dimensions of the observed gene sets.
#' @name helpers
#' @rdname helpers
#' @import asmbPLS
#' @export
#' @returns A list with the observed CIP and GIP distributions.
subsampling_CIP_GIP <- function(object, X.matrix, Y.matrix, K, M, nsubsampling,
                                maxiter, X.dim) {
    table_counts <- table(object@superpathway_input@sample_class)/K
    prop <- as.vector(table_counts[object@superpathway_input@sample_class])
    CIP_GIP_variability <- list()
    for (j in seq_len(nsubsampling)) {
        Nc <- sample(x = min(floor(K/M), 2 * M):K, size = 1)
        # Ensure at least 2 samples per class
        class_indices <- split(
            seq_along(object@superpathway_input@sample_class),
            object@superpathway_input@sample_class)
        min_samples <- unlist(lapply(class_indices,
                                    function(idx) sample(idx, 2)))
        # Determine remaining samples
        remaining_samples <- setdiff(seq_len(K), min_samples)
        Nc_remaining <- Nc - length(min_samples)
        if(Nc_remaining > 0){
            additional_samples <- sample(remaining_samples, size = Nc_remaining,
                                prob = prop[remaining_samples], replace = FALSE)
            training_index <- c(min_samples, additional_samples)
        }else{
            training_index <- min_samples
        }
        # Fit model
        E_matrix_training <- X.matrix[training_index, , drop = FALSE]
        F_matrix_training <- Y.matrix[training_index, , drop = FALSE]
        asmbPLSDA_fit_results <- asmbPLS::asmbPLSDA.fit(
            X.matrix = E_matrix_training, Y.matrix = F_matrix_training,
            PLS.comp = object@hyperparameters_fit@number_PLS, X.dim = X.dim,
            quantile.comb = object@hyperparameters_fit@quantile_comb_table,
            outcome.type = object@hyperparameters_fit@outcome_type,
            center = TRUE, scale = TRUE, maxiter = maxiter)
        object@model_fit$`asmbPLS-DA` <- asmbPLSDA_fit_results
        aux <- CIP_GIP(object)
        # Store results
        if (j == 1) {
            CIP_GIP_variability <- aux
        } else {
            for (b in seq_len(length(X.dim))) {
                CIP_GIP_variability$GIP[[b]] <- cbind(
                    CIP_GIP_variability$GIP[[b]], aux$GIP[[b]])
            }
            CIP_GIP_variability$CIP <- cbind(
                CIP_GIP_variability$CIP, aux$CIP)
        }
    }
    return(CIP_GIP_variability)
}

#' @title Generate Null Distributions by Permutation
#' @description Generate null distributions of CIP and GIP using permutations.
#' @param object Model object.
#' @param X.matrix Predictor matrix.
#' @param Y.matrix Response matrix.
#' @param npermut Number of permutations.
#' @param K Number of samples.
#' @param X.dim Dimensions of the observed gene sets.
#' @param maxiter Maximum iterations.
#' @name helpers
#' @rdname helpers
#' @import asmbPLS
#' @export
#' @returns A list with the null CIP and GIP distributions.
generate_null_distributions <- function(object, X.matrix, Y.matrix,
                                        npermut, K, X.dim, maxiter) {
    NULL_VAR_INF <- list()
    for (j in seq_len((npermut + 1))) {
        X_perm_aux <- permute_X_matrix(X.matrix, K, X.dim)
        # Fit model with permuted data
        Modelpermut <- asmbPLS::asmbPLSDA.fit(
            X.matrix = X_perm_aux, Y.matrix = Y.matrix,
            PLS.comp = object@hyperparameters_fit@number_PLS, X.dim = X.dim,
            quantile.comb = object@hyperparameters_fit@quantile_comb_table,
            outcome.type = object@hyperparameters_fit@outcome_type,
            center = TRUE, scale = TRUE, maxiter = maxiter)
        object@model_fit$`asmbPLS-DA` <- Modelpermut
        aux <- CIP_GIP(object)
        # Store null results
        if (j == 1) {
            NULL_VAR_INF <- aux
        } else {
            for (b in seq_len(length(X.dim))) {
                NULL_VAR_INF$GIP[[b]] <- cbind(
                    NULL_VAR_INF$GIP[[b]], aux$GIP[[b]])
            }
            NULL_VAR_INF$CIP <- cbind(NULL_VAR_INF$CIP, aux$CIP)
        }
    }
    return(NULL_VAR_INF)
}

#' @title Permute X Matrix for Null Distribution
#' @description Permute the X matrix to generate a null distribution.
#' @param X.matrix Predictor matrix.
#' @param K Number of samples.
#' @param X.dim Dimensions of the observed gene sets.
#' @name helpers
#' @rdname helpers
#' @export
#' @returns A permuted X matrix.
permute_X_matrix <- function(X.matrix, K, X.dim) {
    X_perm_aux <- matrix(nrow = K, ncol = ncol(X.matrix))
    colnames(X_perm_aux) <- colnames(X.matrix)
    rownames(X_perm_aux) <- rownames(X.matrix)
    for (b in seq_len(length(X.dim))) {
        if(b == 1){
            indices <- seq(1, X.dim[b])
        }else{
            indices <- (cumsum(X.dim)[b-1]+1):(cumsum(X.dim)[b])
        }
        X_perm_aux[seq_len(K), indices] <- X.matrix[sample(K), sample(indices)]
    }
    return(X_perm_aux)
}

#' @title Calculate P-values for CIP/GIP
#' @description Compute p-values by applying the Mann-Whitney test.
#' @param variability A list of CIP or GIP values for observed distributions.
#' @param null_dist A list of CIP or GIP values for null distributions.
#' @param test_func The test function to use (typically Wilcoxon).
#' @name helpers
#' @rdname helpers
#' @export
#' @returns A data frame of p-values.
calculate_pvalues <- function(variability, null_dist, test_func) {
    lapply(seq_along(variability), function(i) {
        variability_cell <- variability[[i]]
        null_cell <- null_dist[[i]]
        tests <- mapply(function(variability_row, null_row) {
            test_func(variability_row, null_row)
        }, split(variability_cell, row(variability_cell)),
        split(null_cell, row(null_cell)),
        SIMPLIFY = TRUE)
        output <- data.frame("pvalue" = tests)
        rownames(output) <- rownames(variability_cell)
        return(output)
    })
}

#' @title Perform Cross-Validation (LOOCV or K-Fold CV) for asmbPLSDA
#'
#' @description
#' This helper function performs either Leave-One-Out Cross Validation
#' (LOOCV) or K-Fold Cross Validation (KCV) on the given dataset and returns
#' the optimal hyperparameters based on the specified accuracy measure.
#'
#' @param object A `superpathway.input` object containing the data to be used
#' for the cross-validation.
#' @param model_block_matrices A list containing the model block matrices
#' (predictor and response matrices).
#' @param nFC The number of folds for K-fold cross-validation. If `nFC == 1`,
#' LOOCV is performed.
#' @param measure The accuracy measure used to select the optimal model.
#' Default is "F1".
#' @param parallel A logical value indicating whether parallel computation
#' should be used.
#' @param expected_measure_increase Expected decrease in measure per additional
#' PLS component. Default is 0.005.
#' @param maxiter The maximum number of iterations for cross-validation. Default
#' is 100.
#' @param Method The decision rule for prediction (e.g., "fixed_cutoff",
#' "Euclidean_distance_X", etc.). Default is `NULL`.
#' @name helpers
#' @rdname helpers
#' @import asmbPLS
#' @export
#' @returns A list containing the optimal hyperparameters and associated
#' quantile table.
perform_cv <- function(object, model_block_matrices, nFC, measure, parallel,
                        expected_measure_increase, maxiter, Method) {
    nSamples <- as.vector(base::table(object@sample_class))
    if (nFC == 1) {
        message("Running LOOCV")
        return(asmbPLSDA.cv.loo(
            X.matrix = model_block_matrices$block_predictor,
            Y.matrix = model_block_matrices$matrix_response,
            PLS_term = object@hyperparameters_info@number_PLS,
            X.dim = model_block_matrices$block_dim,
            quantile.comb.table =
                object@hyperparameters_info@quantile_comb_table,
            outcome.type = object@hyperparameters_info@outcome_type,
            center = TRUE,
            scale = TRUE,
            measure = measure,
            parallel = parallel,
            expected.measure.increase = expected_measure_increase,
            maxiter = maxiter))
    } else {
        message("Running KCV")
        return(asmbPLS::asmbPLSDA.cv(
            X.matrix = model_block_matrices$block_predictor,
            Y.matrix = model_block_matrices$matrix_response,
            PLS.comp = object@hyperparameters_info@number_PLS,
            X.dim = model_block_matrices$block_dim,
            quantile.comb.table =
                object@hyperparameters_info@quantile_comb_table,
            outcome.type = object@hyperparameters_info@outcome_type,
            k = nFC,
            ncv = 5,
            expected.measure.increase = expected_measure_increase,
            center = TRUE,
            scale = TRUE,
            measure = measure,
            maxiter = maxiter,
            method = Method))
    }
}


#' @title Compute Validation Metrics for the Fitted Model
#'
#' @description
#' This helper function computes various validation metrics, including global
#' significance, CIP/GIP significance, and adjusted p-values for the fitted
#' model based on cross-validation results.
#' @param output The `superpathway.fit.model` object that contains the
#' fitted model and validation information.
#' @param optimal_hyperparameters The optimal hyperparameters obtained
#' from cross-validation.
#' @param model_block_matrices A list containing the model block matrices
#' (predictor and response matrices).
#' @param npermut The number of permutations for significance testing.
#' @param nbObsPermut The number of samples to permute in each permutation.
#' Default is `NULL`.
#' @param maxiter The maximum number of iterations for validation tests.
#' Default is 100.
#' @param global_significance_full Boolean flag indicating whether to return
#' full global significance results.
#' @param CIP.GIP_significance_full Boolean flag indicating whether to return
#' full CIP/GIP significance results.
#' @param type The procedure type for generating CIP/GIP distributions. Can be
#' "jackknife" or "subsampling".
#' @param nsubsampling The number of subsamples for CIP/GIP testing. Default is
#' 100.
#' @param measure The accuracy measure used for validation. Default is "F1".
#' @param Method The decision rule for prediction (e.g., "fixed_cutoff",
#' "Euclidean_distance_X", etc.).
#' @name helpers
#' @rdname helpers
#' @export
#' @returns The updated `superpathway.fit.model` object with the computed
#' validation metrics.
compute_validation_metrics <- function(
        output, optimal_hyperparameters, model_block_matrices, npermut,
        nbObsPermut, maxiter, global_significance_full,
        CIP.GIP_significance_full, type, nsubsampling, measure, Method) {
    Nc <- output@hyperparameters_fit@number_PLS
    CV_error <- optimal_hyperparameters$quantile_table_CV[
        optimal_hyperparameters$optimal_nPLS,
        length(model_block_matrices$block_dim) + get_measure_index(measure)]
    output_global_significance <- permut_asmbplsda(
        output, npermut = npermut, nbObsPermut = nbObsPermut, Nc = Nc, CV_error,
        measure = measure, Method = Method, maxiter = maxiter)
    if (global_significance_full) {
        output@model_validation$`global_significance` <-
            output_global_significance
    } else {
        output@model_validation$`pvalue_global_significance` <-
            output_global_significance$pvalue
    }
    output_CIP_GIP_significance <- CIP_GIP_test(
        output, npermut = npermut, maxiter = maxiter, type = type,
        nsubsampling = nsubsampling)
    if (CIP.GIP_significance_full) {
        output@model_validation$`CIP_significance` <-
            output_CIP_GIP_significance
    } else {
        output@model_validation$`pvalue_CIP_significance` <-
            output_CIP_GIP_significance$`CIP_pvalue`
        output@model_validation$`pvalue_GIP_significance` <-
            output_CIP_GIP_significance$`GIP_pvalue`
    }
    CIP_GIP_adj_pval <- lapply(
        seq_along(output@model_fit$observed_gene_sets),function(j) {
        lambdas <- optimal_hyperparameters$quantile_table_CV[
            seq_len(optimal_hyperparameters$optimal_nPLS),
            seq_len(length(output@model_fit$observed_gene_sets)), drop = FALSE]
        prod_lambdas <- apply(lambdas, 2, function(x) prod(x, na.rm = TRUE))
        m_0 <- ifelse(floor(
            prod_lambdas * lengths(output@model_fit$observed_gene_sets)) == 0,
            1, floor(prod_lambdas *
                        lengths(output@model_fit$observed_gene_sets)))
        pval <- output_CIP_GIP_significance$`GIP_pvalue`[[j]][, 1]
        adj_pval <- pval * m_0[j]
        output_adjpval <- data.frame("adj_p_val" = adj_pval)
        rownames(output_adjpval) <- rownames(
            output_CIP_GIP_significance$`GIP_pvalue`[[j]])
        return(output_adjpval)
    })
    output@model_validation$`adjpvalue_GIP_significance` <- CIP_GIP_adj_pval
    return(output)
}

#' @title Detect gene annotation of a gene set
#' @description
#' For a given gene set it identifies the annotation of the genes, it does so
#' if the genes have more than 50% match with a given annotation. Annotation
#' must be either Ensembl, Entrez or Gene Symbols.
#'
#' @param gene_set A vector with the genes to assess
#' @param mart A `Mart` object from `biomaRt::useMart()`
#'
#' @name helpers
#' @rdname helpers
#' @returns
#' The identified gene annotation or NULL if it was not identified
#' @import checkmate biomaRt
#' @export
#' @examples
#' library(biomaRt)
#' gene_set <- c("IL13", "IL4", "IL5", "IL21")
#' mart <- biomaRt::useMart(biomart = "ensembl",
#' dataset = "hsapiens_gene_ensembl")
#' detect_gene_type(gene_set, mart)
detect_gene_type <- function(gene_set, mart){
    checkmate::assert_class(mart, "Mart")
    possible_filters <- c("ensembl_gene_id", "entrezgene_id",
                            "external_gene_name")
    # Extract organism and if that organism is hsapiens add hgnc_symbol
    # to the possible filters vector
    dataset_name <- base::attributes(mart)$dataset
    organism <- sub("_.*", "", dataset_name)
    if(organism == "hsapiens"){
        possible_filters <- c(possible_filters, "hgnc_symbol")
    }
    matched_filters <- NULL
    for(filter in possible_filters){
        matched_genes <- biomaRt::getBM(attributes = filter, filters = filter,
                                        values = gene_set, mart = mart)
        if(nrow(matched_genes) > length(gene_set) * 0.5){ # At least 50% match
            matched_filters <- filter
            break
        }
    }
    if(is.null(matched_filters)){
        stop("Could not determine gene annotation type. Please provide Ensembl,
                Entrez, or Gene Symbols.")
    }
    return(matched_filters)
}

#' @title Retrieve one to one orthologs
#' @description
#' Retrieves one to one orthologs between from_species and to_species of
#' \link{orthology_mapping}
#'
#' @param annotation A parameter passed from \link{orthology_mapping}, it
#' indicates the annotation of the gene set provided
#' @param gene_set A parameter passed from \link{orthology_mapping}
#' @param mart A parameter passed from \link{orthology_mapping}
#' @param from_species A parameter passed from \link{orthology_mapping}
#' @param to_species A parameter passed from \link{orthology_mapping}
#' @name helpers
#' @rdname helpers
#' @import biomaRt data.table
#' @returns
#' A `data.table` object with the Ensembl identifiers of gene set for
#' from_species and to_species with only one to one orthologs
#' @export
#' @examples
#' annotation <- "external_gene_name"
#' gene_set <- c("IL13", "IL4", "IL5")
#' mart <- biomaRt::useMart(biomart = "ensembl", dataset = paste0("hsapiens",
#' "_gene_ensembl"))
#' retrieve_one2one_orthologs(annotation, gene_set, mart, "hsapiens",
#' "mmusculus")
retrieve_one2one_orthologs <- function(annotation, gene_set, mart,
                                        from_species, to_species){
    # Convert input genes to Ensembl IDs (required for ortholog retrieval)
    gene_conversion <- biomaRt::getBM(
        attributes = c(annotation, "ensembl_gene_id"),
        filters = annotation, values = gene_set,
        mart = mart)
    if(nrow(gene_conversion) == 0) stop("No matching Ensembl IDs found for the
                                        provided genes")
    data.table::setDT(gene_conversion)
    data.table::setnames(gene_conversion, old = annotation, new = "input_gene")
    # Retrieve orthologs using Ensembl IDs
    orthologs <- biomaRt::getBM(
        attributes = c("ensembl_gene_id",
                        paste0(to_species, "_homolog_ensembl_gene"),
                        paste0(to_species, "_homolog_orthology_type")),
        filters = "ensembl_gene_id",
        values = gene_conversion$ensembl_gene_id,
        mart = mart
    )
    data.table::setDT(orthologs)
    data.table::setnames(orthologs,
                        old =c("ensembl_gene_id",
                                paste0(to_species, "_homolog_ensembl_gene"),
                                paste0(to_species, "_homolog_orthology_type")),
                        new = c("from_ensembl", "ortholog", "orthology_type"))
    # Merge back with original input genes
    final_orthologs <- base::merge(gene_conversion, orthologs,
                                    by.x = "ensembl_gene_id",
                                    by.y = "from_ensembl",
                                    all.x = TRUE)
    # Filter for one-to-one orthologs only
    final_orthologs <- final_orthologs[final_orthologs$orthology_type ==
                                        "ortholog_one2one", ]
    return(final_orthologs)
}

#' @title Translate Fold Changes onto human gene expression
#' @description
#' Applies the biological link function conditions onto a predictor block
#' matrix. The resulting gene expression of the predictor block are the
#' cases defined in the biological link function.
#' @param model_object A \link{superpathway.fit.model-class} object
#' @param b A parameter passed from \link{singIST_treat}. The index of
#' current iteration block.
#' @param samples A parameter passed from \link{singIST_treat}. The samples
#' to modify its gene expression from `predictor_block`
#' @param predictor_block A parameter passed from \link{singIST_treat}.
#' The predictor block of matrices from asmbPLSDA to modify its gene expression.
#' @param FC A parameter passed from \link{singIST_treat}. A `data.frame` with
#' the Fold Changes, for a cell type, of each gene.
#' @name helpers
#' @rdname helpers
#' @returns
#' The predictor block matrix updated with the FC translation
#' @export
FCtoExpression <- function(model_object, b, samples, predictor_block, FC){
    indices <- which(colnames(predictor_block) %in% rownames(FC))
    indices_FC <- match(colnames(predictor_block[, indices, drop = FALSE]),
                        rownames(FC))
    mu <- model_object@model_fit$`asmbPLS-DA`$X_col_mean[,indices, drop = FALSE]
    r <- t(FC[indices_FC, "avg_log2FC", drop = FALSE])
    mu_per_r <- mu*r
    min_col <- apply(predictor_block[samples, indices, drop = FALSE], 2, min)
    # Check condition C = min{x + r*mu} >= 0
    sum <- min_col + mu_per_r
    check_negativity_condition <-  apply(sum, 2,
                                            FUN = function(x)
                                            return(ifelse(x >= 0, TRUE, FALSE))
    )
    # Translate FC onto gene expression
    for(i in seq(1, ncol(mu))){
        if(check_negativity_condition[i]){
            predictor_block[samples, indices[i]] <-
                predictor_block[samples, indices[i]] + mu_per_r[1, i]
        }else{
            predictor_block[samples, indices[i]] <-
                predictor_block[samples, indices[i]] - min_col[i]
        }
    }
    return(predictor_block)
}

#' @title Center and scale predictor block matrices
#' @description
#' Centers and scales each column of the predictor block matrices. The
#' centering and scaling is according to the centroid and variance estimated in
#' `fit_asmb`.
#' @param data Block of predictor matrices to center and scale
#' @param fit_asmb An asmbPLSDA fitted model
#' @name helpers
#' @rdname helpers
#' @returns
#' The object `data` centered and scaled.
#' @export
center_scale <- function(data, fit_asmb){
    center <- fit_asmb$center
    scale <- fit_asmb$scale
    X.dim <- fit_asmb$X_dim
    X_col_mean <- fit_asmb$X_col_mean
    X_col_sd <- fit_asmb$X_col_sd
    if(center){
        if(scale){
            for(i in seq(1,sum(X.dim))){
                data[, i] <- (data[, i] - X_col_mean[i])/X_col_sd[i]
            }
        }else{
            for(i in seq(1,sum(X.dim))){
                data[, i] <- data[, i] - X_col_mean[i]
            }
        }
    }else{
        if(scale){
            for(i in seq(1,sum(X.dim))){
                data[, i] <- data[i, ]/X_col_sd[i]
            }
        }
    }
    return(data)
}

#' @title Extract indices for a block
#' @description
#' Given a block and the dimensions of all blocks it returns the indices of the
#' genes belonging to that block within the predictor block matrix
#'
#' @param j Block to return indices for
#' @param X.dim Vector with number of genes of each block
#'
#' @name helpers
#' @rdname helpers
#'
#' @returns
#' A vector with the indices of the predictor block matrix for the requestes
#' block
#' @export
#' @examples
#' X.dim <- c(30,40,60)
#' j <- 2
#' get_indices(j, X.dim)
get_indices <- function(j, X.dim){
    if(j == 1){
        indices <- seq(1, X.dim[j])
    }else{
        indices <- (cumsum(X.dim)[j-1]+1):(cumsum(X.dim)[j])
    }
    return(indices)
}

#' @title Deflate prediction
#' @description
#' Performs loading deflation for a given predictor block and PLS component
#'
#' @param data Matrix of predictor block to deflate
#' @param PLS Numeric value indicating the PLS component
#' @param delta_cbind Gene contributions (loadings) used to deflate the blocks
#' @param fit_asmb asmbPLSDA fitted model
#' @name helpers
#' @rdname helpers
#' @returns
#' The `data` matrix loading deflated
#' @export
deflate_prediction <- function(data, PLS, delta_cbind, fit_asmb){
    X.dim <- fit_asmb$X_dim
    for(j in seq(1,length(X.dim))){
        indices <- get_indices(j, X.dim)
        C_b <- as.matrix(data[, indices])
        # Deflation by block loadings
        p_temp <- as.matrix(fit_asmb$X_loading[[j]][,PLS])
        p <- as.matrix(delta_cbind[[j]][,PLS])
        data[, indices] <- C_b - p*p_temp
    }
    return(data)
}
