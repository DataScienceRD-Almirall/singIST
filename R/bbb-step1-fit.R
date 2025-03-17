#' @title Build predictor and response blocks with superpathway input
#'
#' @description
#' Builds the predictor block and matrix response for its fit in asmbPLS-DA
#'
#' @param object \link{superpathway.input-class} object
#' @rdname matrixToBlock-method
#'
#' @returns A list containing the predictor block, response matrix, dimension of
#' each block and observed gene sets with respect to gene_sets_celltype for the
#' original pseudobulk_lognorm matrix, for its use in asmbPLS-DA fit
#'
#' @import checkmate stats
#' @export
#' @examples
#' data <- example_superpathway_input
#' matrixToBlock(data)
matrixToBlock.superpathway.input <- function(object){
    checkmate::assert_class(object, "superpathway.input")
    matrix <- object@pseudobulk_lognorm
    gene_sets_celltype <- object@superpathway_info@gene_sets_celltype
    observed_gene_sets <- base::lapply(gene_sets_celltype, function(x)
                                        intersect(x, colnames(matrix)))
    block_dim <- vapply(observed_gene_sets, length, FUN.VALUE = numeric(1))
    if(all(block_dim > 0)){
        split <- base::do.call(base::rbind, base::strsplit(rownames(matrix),
                                                            "_"))
        block_celltypes <- split[, 1]
        block_sample_id <- split[, 2]
        block_predictor <- matrix(ncol = sum(block_dim),
                                    nrow = length(unique(block_sample_id)))
        rownames(block_predictor) <- unique(block_sample_id)
        colnames(block_predictor) <- as.vector(unlist(Map(
            function(prefix, value)
            paste0(prefix, "*", value),unique(block_celltypes),
            observed_gene_sets)))
    # We use Reduce() to update block_predictor accumulatively
    block_predictor <- base::Reduce(function(bm, i){
        update_block(unique(block_celltypes)[i], observed_gene_sets[[i]],
                        bm, matrix = matrix)
        }, seq_along(unique(block_celltypes)), init = block_predictor)
    }else{
        stop("There is at least one cell type gene set that is void check that
                at least one gene exists in your pseudobulk_lognorm matrix for
                all cell type gene sets.")
        }
    # Build predictor block
    categories_class <- base::factor(object@sample_class)
    categories_class <- stats::relevel(categories_class,ref = object@base_class)
    matrix_response <- stats::model.matrix(~ categories_class -1)
    if(object@hyperparameters_info@outcome_type == "binary"){
        # Remove reference class column
        matrix_response <- matrix_response[ , 2:ncol(matrix_response),
                                            drop = FALSE]
        }else{
            matrix_response <- matrix_response[ , , drop = FALSE]
            }
    # Return predictor_block, matrix_response and block_dim
    output <- list(block_predictor = block_predictor,
                    matrix_response = matrix_response,block_dim = block_dim,
                    observed_gene_sets = observed_gene_sets)
    return(output)
}

#' @title Compute performance metrics of predicted asmbPLSDA
#'
#' @param Y_predict Predicted matrix from asmbPLSDA
#' @param Y_true True class used to fit asmbPLSDA
#' @param outcome.type Outcome type either `"binary"` or `"multiclass"`
#'
#' @returns
#' A vector with accuracy, balanced accuracy, precision, recall and F1 metric
#' @export
#' @examples
#' Results_comparison_measure(c(1,0,1,0,1), c(0,0,1,1,1),
#' outcome.type = "binary")
Results_comparison_measure <- function(Y_predict,
                                        Y_true,
                                        outcome.type= c("binary","multiclass")){
    Y_col <- length(as.vector(Y_true))
    n_match <- n_TP <- n_TN <- n_FP <- n_FN <-
    balanced_accuracy_multicalss <- n_recall_multiclass <- 0
    for (i in seq_len(Y_col)) {
        temp_accu <- which(as.vector(Y_predict)[i] == as.vector(Y_true)[i])
        temp_TP <- which(as.vector(Y_predict)[i] == 1 &&
                            as.vector(Y_true)[i] == 1)
        temp_TN <- which(as.vector(Y_predict)[i] == 0 &&
                            as.vector(Y_true)[i] == 0)
        temp_FP <- which(as.vector(Y_predict)[i] == 1 &&
                            as.vector(Y_true)[i] == 0)
        temp_FN <- which(as.vector(Y_predict)[i] == 0 &&
                            as.vector(Y_true)[i] == 1)
        n_temp_TP <- length(temp_TP)
        n_temp_TN <- length(temp_TN)
        n_temp_FP <- length(temp_FP)
        n_temp_FN <- length(temp_FN)
        n_match <- n_match + temp_accu
        n_TP <- n_TP + n_temp_TP
        n_TN <- n_TN + n_temp_TN
        n_FP <- n_FP + n_temp_FP
        n_FN <- n_FN + n_temp_FN
        n_recall_multiclass <- n_temp_TP/(n_temp_TP + n_temp_FN)
        balanced_accuracy_multicalss <- balanced_accuracy_multicalss +
            n_recall_multiclass
    }
    accuracy <- (n_TP + n_TN)/(n_TP + n_TN + n_FP + n_FN)
    precision <- n_TP/(n_TP + n_FP)
    recall <- n_TP/(n_TP + n_FN) # Sensitivity also
    specificity <- n_TN/(n_TN + n_FP)
    F1 <- 2 * precision * recall/(precision + recall)
    balanced_accuracy <- 0
    balanced_accuracy_binary <- (recall + specificity)/2
    balanced_accuracy_multicalss <- balanced_accuracy_multicalss/Y_col
    if (outcome.type == "binary") {
        balanced_accuracy <- balanced_accuracy_binary
    }else if(outcome.type == "multiclass"){
        balanced_accuracy <- balanced_accuracy_multicalss
    }
    output <- c("accuracy" = accuracy, "balanced_accuracy" = balanced_accuracy,
                "precision" = precision, "recall" = recall, "F1" = F1)
    return(output)
}

#' @title Leave-one-out Cross-validation
#'
#' @param X.matrix Predictor block matrix from \code{matrixToBlock}
#' @param Y.matrix Response matrix from \code{matrixToBlock}
#' @param PLS_term An integer with the number of PLS components to use passed
#' from \link{hyperparameters-class} obect
#' @param X.dim A list with the observed gene set size for each cell type
#' from \code{matrixToBlock}
#' @param quantile.comb.table A matrix with the quantile comb table passed
#' from \link{hyperparameters-class} object
#' @param outcome.type A character indicating `binary` or `multiclass` passed
#' from \link{hyperparameters-class} object
#' @param Method A parameter passed from \code{fitOptimal}
#' @param measure A parameter passed from \code{fitOptimal}
#' @param parallel A parameter passed from \code{fitOptimal}
#' @param cores A parameter passed from \code{fitOptimal}
#' @param expected.measure.increase A parameter passed from \code{fitOptimal}
#' @param center A parameter passed from \code{fitOptimal}
#' @param scale A parameter passed from \code{fitOptimal}
#' @param maxiter A parameter passed from \code{fitOptimal}
#'
#' @returns
#' A list containing the optimal quantiles for each PLS component and the
#' optimal number of PLS components.
#' @export
#'
#' @examples
#' data <- example_superpathway_input
#' matrices <- matrixToBlock(data)
#' X.matrix <- matrices$block_predictor
#' Y.matrix <- matrices$matrix_response
#' X.dim <- matrices$block_dim
#' quantile.comb.table <- slot(slot(data, "hyperparameters_info"),
#' "quantile_comb_table")
#' outcome.type <- slot(slot(data, "hyperparameters_info"), "outcome_type")
#' asmbPLSDA.cv.loo(X.matrix, Y.matrix, PLS_term = 1, X.dim,quantile.comb.table,
#' Method = NULL, measure = "B_accuracy", parallel = TRUE, cores = NULL,
#' outcome.type = outcome.type, expected.measure.increase = 0.005,
#' center = TRUE, scale = TRUE,maxiter = 100)
asmbPLSDA.cv.loo <- function(X.matrix, Y.matrix, PLS_term = 1, X.dim,
                                quantile.comb.table, outcome.type =
                                c("binary", "multiclass"), Method = NULL,
                                measure = "B_accuracy", parallel = FALSE,
                                cores = NULL, expected.measure.increase = 0.005,
                                center = TRUE, scale = TRUE, maxiter = 100){
    n_group <- ncol(Y.matrix)
    measure_selected <- get_measure_index(measure)
    K <- nrow(Y.matrix)
    n_quantile_comb <- nrow(quantile.comb.table)
    quantile_table_CV <- matrix(data = 0, nrow = PLS_term,
                                ncol = (length(X.dim) + 5))
    for (i in seq_len(PLS_term)) {
        results_CV_summary_n <- matrix(0, nrow = n_quantile_comb, ncol = K)
        F_matrix_validation_bind <- matrix(0, nrow = n_quantile_comb, ncol = K)
        if (parallel) {
            results <- execute_parallel_cv(
                K, cores, results_CV_summary_n, F_matrix_validation_bind,
                X.matrix, Y.matrix, i, X.dim, quantile.comb.table, outcome.type,
                quantile_table_CV, Method, measure, expected.measure.increase,
                center, scale, maxiter)
            results_CV_summary_n <- results$results_CV_summary_n
            F_matrix_validation_bind <- results$F_matrix_validation_bind
        } else {
            results <- execute_sequential_cv(
                K, n_quantile_comb, results_CV_summary_n,
                F_matrix_validation_bind,X.matrix, Y.matrix, i, X.dim,
                quantile.comb.table, outcome.type, quantile_table_CV,measure,
                expected.measure.increase, center, scale, maxiter, Method)
            results_CV_summary_n <- results$results_CV_summary_n
            F_matrix_validation_bind <- results$F_matrix_validation_bind
        }
        measure_acc <- performance_measures(n_quantile_comb,
                                            results_CV_summary_n,
                                            F_matrix_validation_bind,
                                            outcome.type,
                                            measure_selected)
        index_max_measure <- which.max(measure_acc)
        quantile_table_CV[i, seq_len(length(X.dim))] <-
            quantile.comb.table[index_max_measure, ]
        quantile_table_CV <- compute_final_measures(
            K, X.matrix, Y.matrix, i, X.dim, quantile_table_CV, outcome.type,
            center, scale, maxiter, Method)
    }
    optimal_nPLS <- select_optimal_PLS(PLS_term, quantile_table_CV, X.dim,
                                    measure_selected, expected.measure.increase)
    return(list("quantile_table_CV" = quantile_table_CV,
                "optimal_nPLS" = optimal_nPLS))
}

#' @title Compute Cell Importance Projection (CIP) and Gene Importance
#' Projection (GIP)
#' @description
#' Computes CIP and GIP metrics from a \link{superpathway.fit.model-class}
#' object for the target class
#'
#' @param object A \link{superpathway.fit.model-class} object
#'
#' @returns
#' A list with the CIP and GIP metrics for all cell types. The metrics are
#' computed for the target class.
#' @import checkmate stringr
#' @export
#' @examples
#' CIP_GIP(example_superpathway_fit_model)
CIP_GIP <- function(object){
    checkmate::assert_class(object, "superpathway.fit.model")
    # Identify target class position
    class_position <- stringr::str_replace(
        colnames(object@model_fit$response_matrix), ".*categories_class", "")
    target_class_position <- which(
        class_position == object@superpathway_input@target_class)
    # asmbPLSDA model
    model <- object@model_fit$`asmbPLS-DA`
    # Loadings
    w_super <- model$X_super_weight^2
    q <- model$Y_weight[target_class_position, , drop = FALSE]
    nblocks <- length(object@model_fit$observed_gene_sets)
    # Compute CIP
    CIP <- (w_super %*% t(q))/sum(q)
    rownames(CIP) <- object@superpathway_input@superpathway_info@celltypes
    # Compute GIP
    GIP <- vector("list", nblocks)
    names(GIP) <- object@superpathway_input@superpathway_info@celltypes
    for(b in seq_len(nblocks)){
        w <- model$X_weight[[b]]^2
        GIP[[b]] <- (w %*% t(q))/sum(q)
        rownames(GIP[[b]]) <- rownames(w)
    }
    return(list("GIP" = GIP, "CIP" = CIP))
}

#' @title Mann-Whitney Wilcoxon test p-value
#'
#' @param ref_distr A vector with the reference distribution
#' @param null_distr A vector with the null distribution
#'
#' @returns
#' A pvalue with the Mann-Whitney Wilcoxon test with the "greater" as the
#' alternative hypothesis
#'
#' @import stats
#' @export
#' @examples
#' ref_distr <- rnorm(100, mean = 30, sd = 2)
#' null_distr <- rnorm(100, mean = 0, sd = 1)
#' wilcox_CIP_GIP(ref_distr, null_distr)
wilcox_CIP_GIP <- function(ref_distr, null_distr){
    return(stats::wilcox.test(
        ref_distr, null_distr, alternative = "greater", exact = FALSE)$p.value)
}

#' @title Permutation test for asmbPLSDA global significance
#'
#' @description
#' Performs permutation testing for asmbPLS-DA to evaluate model validity.
#'
#' @param object A \link{superpathway.fit.model-class} object.
#' @param npermut Number of permutations (default: 100).
#' @param nbObsPermut Number of samples to permute per iteration
#' (default: NULL).
#' @param Nc Number of samples dropped per permutation (default: 1).
#' @param CV_error Cross-validation error of the fitted model.
#' @param measure Accuracy measure (`"F1"`, `"accuracy"`, `"B_accuracy"`,
#' `"precision"`, `"recall"`, default: `"B_accuracy"`).
#' @param Method Decision rule for prediction (default: NULL).
#' @param maxiter Maximum iterations (default: 100).
#'
#' @return A list with permutation statistics, p-value, and confidence
#' intervals.
#' @export
#' @examples
#' permut_asmbplsda(example_superpathway_fit_model, npermut = 5, Nc = 1,
#' CV_error = 1)
permut_asmbplsda <- function(object, npermut = 100, nbObsPermut = NULL,
                                Nc = 1, CV_error, measure = "B_accuracy",
                                Method = NULL, maxiter = 100) {
    Y.matrix <- object@model_fit$`asmbPLS-DA`$Y_group
    X.matrix <- object@model_fit$predictor_block
    nr <- nrow(Y.matrix)
    q <- ncol(Y.matrix)
    res <- initialize_results(npermut, q)
    for (j in seq_len(npermut + 1)) {
        Ypermut <- permute_Y_matrix(Y.matrix, nr, nbObsPermut, j)
        res <- compute_permutation_stats(res, Y.matrix, Ypermut, j, q, nr)
        s <- select_samples(object, nr, Nc)
        X_train <- X.matrix[-s, , drop = FALSE]
        X_val <- X.matrix[s, , drop = FALSE]
        Y_train <- Ypermut[-s, , drop = FALSE]
        Y_val <- Ypermut[s, , drop = FALSE]
        Modelpermut <- fit_permuted_model(object, X_train, Y_train, maxiter)
        res <- evaluate_performance(res, Modelpermut, X_train, X_val,
                                    Y.matrix, s, measure, j, nr, Method, object)
    }
    null_errors <- as.vector(
        res$prct.Ychange.values[ , ncol(res$prct.Ychange.values)])
    # If F1 is NaN => Recall + Precision = 0, the model performance is poor
    # impute to 0
    if(measure == "F1"){
        null_errors[is.nan(null_errors)] <- 0
    }
    res$pvalue <- compute_pvalue(null_errors, CV_error)
    res$IC <- compute_IC95(null_errors)
    return(res)
}


#' @title Cell and Gene Importance Projections statistical significance
#'
#' @description
#' Computes Cell and Gene Importance Projection observed distribution from
#' fitted asmbPLSDA, and its associated null distributions by permuting the
#' block of predictor matrices. Returns a pvalue of the Mann-Whitney Wilcoxon
#' between the observed and null distribution for each CIP and GIP.
#'
#' @param object A \link{superpathway.fit.model-class} object
#' @param npermut Number of permutations on response block matrices
#' @param maxiter An integer indicating the maximum number of iterations.
#' If `NULL` the default is 100.
#' @param type Either `jackknife` or `subsampling`. If `jackknife` then the CIP
#' and GIP observed distribution is generated by a jackknife procedure. If
#' `subsampling` the CIP and GIP observed distribution is generated by
#' subsampling the number of samples without replacement, each subsample is
#' guaranteed to contain at least 2 samples per class. If a LOOCV was performed
#' or one has small sample size it is recommended to select `jackknife`,
#' otherwise select `subsampling`.
#' @param nsubsampling Number of subsamples to generate CIP and GIP observed
#' distributions. By default 100.
#'
#' @import asmbPLS checkmate
#' @returns
#' A list containing; observed distributions of CIP and GIP (variability_param);
#' its associated null distributions generated by permutations (NULL_CIP_GIP);
#' the unadjusted pvalue of Mann-Whitney Wilcoxon for CIP distribution
#' (CIP_pvalue); and for GIP distribution (GIP_pvalue).
#' @export
#' @examples
#' CIP_GIP_test(example_superpathway_fit_model, npermut = 3, type = "jackknife")
CIP_GIP_test <- function(object, npermut = 100, maxiter = 100,
                        type = c("jackknife", "subsampling"),
                        nsubsampling = 100) {
    checkmate::assert_choice(type, choices = c("jackknife", "subsampling"))
    # Extract data from the object
    K <- nrow(object@model_fit$`asmbPLS-DA`$Y_group)
    M <- length(unique(object@superpathway_input@sample_class))
    X.matrix <- object@model_fit$predictor_block
    Y.matrix <- object@model_fit$response_matrix
    X.dim <- lengths(object@model_fit$observed_gene_sets)
    # Initialize result containers
    CIP_GIP_variability <- NULL_VAR_INF <- list()
    # Calculate CIP/GIP distributions based on resampling type
    if (type == "jackknife") {
        CIP_GIP_variability <- jackknife_CIP_GIP(object, X.matrix, Y.matrix, K,
                                                    maxiter, X.dim)
    } else {
        CIP_GIP_variability <- subsampling_CIP_GIP(object, X.matrix, Y.matrix,
                                                    K, M, nsubsampling, maxiter,
                                                    X.dim)
    }
    # Generate null distributions by permutation
    NULL_VAR_INF <- generate_null_distributions(object, X.matrix, Y.matrix,
                                                npermut, K, X.dim, maxiter)
    # Compute p-values for CIP and GIP distributions
    GIP_pvalue <- calculate_pvalues(CIP_GIP_variability$GIP, NULL_VAR_INF$GIP,
                                    wilcox_CIP_GIP)
    CIP_pvalue <- lapply(seq_along(CIP_GIP_variability$GIP), function(i){
        variability_cell_CIP <- CIP_GIP_variability$CIP[i,]
        null_cell_CIP <- NULL_VAR_INF$CIP[i,]
        tests <- wilcox_CIP_GIP(variability_cell_CIP, null_cell_CIP)
        output <- tests
        return(output)
    })
    # Return results
    return(list("variability_param" = CIP_GIP_variability,
                "NULL_CIP_GIP" = NULL_VAR_INF, "CIP_pvalue" = CIP_pvalue,
                "GIP_pvalue" = GIP_pvalue))
}

#' @title Cross validation and fit of asmbPLSDA
#'
#' @description
#' Performs Cross Validation of the provided superpathway input, fits the
#' optimal model and computes its validation metrics. The Cross Validation can
#' either be Leave One-Out Cross Validation (LOOCV) or K-Fold Cross Validation
#' (KCV). A LOOCV is performed if the number of folds was set to 1 or if the
#' number of samples per class is less than 3 for any class. A K-Fold Cross
#' Validation (KCV) is performed if the number of folds is greater or equal
#' than 3 and the number of samples per class is always greater than the number
#' of folds. If the number of samples is low for some of the classes LOOCV is
#' recommended.
#'
#' @param object A \link{superpathway.input-class} object to fit optimal
#' asmbPLSDA.
#' @param parallel A boolean indicating whether to parallelize (`TRUE`)
#' for LOOCV on quantile combination or not (`FALSE`). Note this option is only
#' available for LOOCV and not KCV. Default is `FALSE`.
#' @param cores Integer number of cores to use if parallize is `TRUE`. By
#' default if `NULL` cores are assigned to
#' `parallel::detectCores(logical=FALSE)-1`.
#' @param measure Accuracy measure to be used to select optimal asmbPLSDA model.
#' Default is F1 measure. Options are: F1, accuracy, B_accuracy, precision
#' and recall.
#' @param Method Decision rule used for prediction. For binary outcome
#' `fixed_cutoff` (default), `Euclidean_distance_X`, and
#' `Mahalanobis_distance_X`. For categorical otcome with more than 2 levels,
#' the methods include `Max_Y` (default), `Euclidean_distance_X`,
#' `Mahalanobis_distance_X`, `Euclidean_distance_Y`, and
#' `PCA_Mahalanobis_distance_Y`. If `NULL` the default method is used for the
#' respective outcome binary.
#' @param expected_measure_increase A double indicating the measure you expect
#' to decrease by percent after including one more PLS component, this will
#' affect the selection of optimal number of PLS components. If `NULL` the
#' default is 0.005 (0.5%).
#' @param maxiter An integer indicating the maximum number of iterations.
#' If `NULL` the default is 100.
#' @param global_significance_full A boolean indicating whether to return a list
#' with information of each permutation for the global
#' significance test of asmbPLSDA. By default `FALSE`. Note that if the number
#' of permutations that is set is large, storing this information can
#' be a burden on memory.
#' @param CIP.GIP_significance_full A boolean indicating whether to return a
#' list with the observed and null distributions of CIP and GIP or only the
#' pvalue and adjusted pvalue. By default `FALSE`. Note that if the number of
#' permutations that is set is large, storing this information can be a burden
#' on memory.
#' @param npermut Number of permutations for the tests. By default 100.
#' Parameter passed onto \link{permut_asmbplsda} and \link{CIP_GIP_test}.
#' @param nbObsPermut An integer indicating the number of samples to permute
#' in each permutation. By default `NULL`. If `NULL` the number of samples to
#' permute at each permutation is randomly chosen (for each permutation).
#' Parameter passed onto \link{permut_asmbplsda}.
#' @param type Either `jackknife` or `subsampling`. If `jackknife` then the CIP
#' and GIP observed distribution is
#' generated by a jackknife procedure. If `subsampling` the CIP and GIP observed
#' distribution is generated by subsampling the number of samples without
#' replacement, each subsample is guaranteed to contain at least 2 samples per
#' class. If a LOOCV was performed or one has small sample size it is
#' recommended to select `jackknife`, otherwise select `subsampling`.
#' Passed onto \link{CIP_GIP_test}.
#' @param nsubsampling Number of subsamples to generate CIP and GIP observed
#' distributions. By default 100. Passed onto \link{CIP_GIP_test}.
#' @import asmbPLS checkmate
#' @rdname fitOptimal-method
#'
#' @returns
#' A \link{superpathway.fit.model-class} object with; a
#' \link{superpathway.input-class} object used for CV and model fit;
#' a \link{hyperparameters-class} object with the hyperparameters used to fit
#' the optimal model (includes optimal quantiles and PLS components from the
#' CV step); a list with the fitted model information including: predictor
#' and response matrices, observed gene sets, from \code{matrixToBlock},
#' and asmbPLSDA output; a list with the validaton metrics of fitted model.
#' @examples
#' # fitOptimal with jackknife for CIP/GIP statistics and 10 permutations
#' # for the global significance test of the optimal model
#' fitOptimal(example_superpathway_input, npermut = 10, type = "jackknife")
#' # fitOptimal with subsampling for CIP/GIP statistics with
#' # 10 subsamples and 50 permutations for the global significance test of the
#' # optimal model
#' fitOptimal(example_superpathway_input, npermut = 50, type = "subsampling",
#' nsubsampling = 10)
fitOptimal.superpathway.input <- function(
        object, parallel = FALSE, cores = NULL, measure = "B_accuracy",
        Method = NULL, expected_measure_increase = 0.005, maxiter = 100,
        global_significance_full = FALSE, CIP.GIP_significance_full = FALSE,
        npermut = 100, nbObsPermut = NULL, type = "jackknife",
        nsubsampling = 100) {
    output <- new("superpathway.fit.model", superpathway_input = object,
                    hyperparameters_fit = object@hyperparameters_info,
                    model_fit = list(), model_validation = list())
    measure_selected <- get_measure_index(measure)
    if (min(as.vector(base::table(object@sample_class))) <= 1) {
        stop("At least one class has 1 or fewer samples.")
    }
    model_block_matrices <- matrixToBlock(object)
    nFC <- ifelse(is.null(object@hyperparameters_info@folds_CV), 5,
                    object@hyperparameters_info@folds_CV)
    if(!(min(as.vector(base::table(object@sample_class))) >= nFC)){
        if(!(min(as.vector(base::table(object@sample_class)))>= as.integer(3))){
            message("Cannot run KCV with initial Folds, running LOOCV instead")
            nFC <- as.integer(1)
        }else{
            message("Cannot run KCV with initial Folds, run KCV with 3-Folds")
            nFC <- as.integer(3)
        }
    }
    optimal_hyperparameters <- perform_cv(object, model_block_matrices, nFC,
        measure, parallel, cores, expected_measure_increase, maxiter, Method)
    output@hyperparameters_fit@number_PLS <- as.integer(
        optimal_hyperparameters$optimal_nPLS)
    output@hyperparameters_fit@quantile_comb_table <-
        optimal_hyperparameters$quantile_table_CV
    output@hyperparameters_fit@folds_CV <- as.integer(nFC)
    optimal_fit <- asmbPLS::asmbPLSDA.fit(
        X.matrix = model_block_matrices$block_predictor,
        Y.matrix = model_block_matrices$matrix_response,
        PLS.comp = output@hyperparameters_fit@number_PLS,
        X.dim = model_block_matrices$block_dim, center = TRUE, scale = TRUE,
        quantile.comb = output@hyperparameters_fit@quantile_comb_table,
        outcome.type = output@hyperparameters_fit@outcome_type)
    output@model_fit <- list(
        "predictor_block" = model_block_matrices$block_predictor,
        "response_matrix" = model_block_matrices$matrix_response,
        "observed_gene_sets" = model_block_matrices$observed_gene_sets,
        "asmbPLS-DA" = optimal_fit)
    output <- compute_validation_metrics(output, optimal_hyperparameters,
            model_block_matrices, npermut, nbObsPermut, maxiter,
            global_significance_full, CIP.GIP_significance_full, type,
            nsubsampling, measure, Method)
    return(output)
}
