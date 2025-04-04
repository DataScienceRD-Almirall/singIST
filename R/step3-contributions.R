#' @title Compute predictor scores
#' @description
#' Computes scores from the predictor to later derive the superpathway's score,
#' cell type contribution and gene contributions, for the target class
#'
#' @param object A \link{superpathway.fit.model-class} object passed from
#' \link{derive_contributions}
#' @param data Block of predictor matrices to compute scores from
#' @param sample Current sample from `data` to compute scores
#'
#' @returns
#' A list containing the needed parameters to compute superpathway's score,
#' cell type contributions and needed scores to compute gene contributions
#' @export
#' @examples
#' data(example_mapping_organism)
#' mapped <- example_mapping_organism
#' data(example_superpathway_fit_model)
#' model <- example_superpathway_fit_model
#' singIST_samples <- biological_link_function(mapped, model)$singIST_samples
#' # Derive the scores for sample 2
#' derive_scores(model, singIST_samples, 2)
derive_scores <- function(object, data, sample){
    fit_asmb <- object@model_fit$`asmbPLS-DA`
    delta_cbind <- Delta <- gamma <- Gamma <- c()
    Y_fit <- 0
    n.PLS <- ncol(fit_asmb$Y_weight)
    FC_applied_test <- as.matrix(t(data[sample, colnames(data)]))
    X.dim <- fit_asmb$X_dim
    Y_col_mean <- fit_asmb$Y_col_mean
    X_col_mean <- fit_asmb$X_col_mean
    X_col_sd <- fit_asmb$X_col_sd
    FC_applied_test <- center_scale(FC_applied_test, fit_asmb)
    cell_types <- object@superpathway_input@superpathway_info@celltypes
    class_position <- stringr::str_replace(
        colnames(object@model_fit$response_matrix), ".*categories_class", "")
    target_class_position <- which(
        class_position == object@superpathway_input@target_class)
    for(i in seq(1, n.PLS)){
        for(j in seq(1, length(X.dim))){
            indices <- get_indices(j, X.dim)
            C_b <- as.matrix(FC_applied_test[, indices])
            value_aux <- unique(sub("\\*.*", "", rownames(C_b)))
            index_aux <- which(cell_types == value_aux)
            omega <- as.matrix(fit_asmb$X_weight[[j]][,i])
            delta <- C_b*omega/(sqrt(X.dim[j]))
            rownames(delta) <- sub(".*\\*", "", rownames(delta))
            colnames(delta) <- paste0(cell_types[index_aux], "_",i)
            if(i == 1){delta_cbind[[j]] <- delta}
            if(i > 1){delta_cbind[[j]] <- cbind(delta_cbind[[j]], delta)}
        }
        Delta <- as.matrix(cbind(Delta, stats::setNames(vapply(delta_cbind,
            function(x)colSums(x[, i, drop = FALSE]), FUN.VALUE = 1), vapply(
            delta_cbind, function(x)colnames(x)[i], FUN.VALUE = "character"))))
        rownames(Delta) <- sub("\\_.*", "", rownames(Delta))
        omega_super <- as.matrix(fit_asmb$X_super_weight[,i])
        aux <- cbind(Delta[,i], omega_super)
        gamma <- cbind(gamma, aux[,1]*aux[,2])
        Gamma <- cbind(Gamma, sum(gamma[,i]))
        Gamma <- as.matrix(Gamma)
        FC_applied_test <- deflate_prediction(FC_applied_test, i,
                                                delta_cbind,fit_asmb)
        q <- fit_asmb$Y_weight[target_class_position, i]
        Y_fit <- Y_fit + Gamma[,i]*q
    }
    if(fit_asmb$center){Y_fit <- Y_fit + Y_col_mean[target_class_position]}
    output <- list(
        Delta = Delta, delta = delta_cbind, gamma = gamma, Gamma = Gamma,
        Y_pred_num = Y_fit, Y_weight = fit_asmb$Y_weight[target_class_position,
                                                            , drop = FALSE])
    return(output)
}

#' @title Derive superpathway score, cell type contribution and gene
#' contribution
#' @description
#' Computes the superpathway score, its cell type contribution and gene
#' contribution for a block of predictor matrices for its later use to compute
#' recapitulations
#' @param model_object A \link{superpathway.fit.model-class} object with the
#' fitted asmbPLSDA
#' @param data A matrix with the block of predictor matrices to compute
#' score and contributions from
#' @import checkmate
#' @returns
#' A list with the superpathway score for each sample, cell type and gene
#' contributions to the former
#' @export
#' @examples
#' data(example_mapping_organism)
#' mapped <- example_mapping_organism
#' data(example_superpathway_fit_model)
#' model <- example_superpathway_fit_model
#' singIST_samples <- biological_link_function(mapped, model)$singIST_samples
#' derive_contributions(model, singIST_samples)
derive_contributions <- function(model_object, data){
    checkmate::assert_class(model_object, "superpathway.fit.model")
    superpathway_score <- celltype_contribution <- gene_contribution <- c()
    fit_asmb <- model_object@model_fit$`asmbPLS-DA`
    # Derive scores for each sample under analysis
    for(i in seq(1, nrow(data))){
        output <- derive_scores(model_object, data, i)
        # Superpathway score
        superpathway_score <- cbind(superpathway_score, output$Y_pred_num)
        # Cell type contribution to superpathway score
        celltype_contribution <- cbind(celltype_contribution, output$gamma %*%
                                            as.matrix(t(output$Y_weight)))
        # Gene contribution to cell type recapitulation
        j_num <- 1
        for(j in rownames(output$gamma)){
            w_super <- as.matrix(fit_asmb$X_super_weight[j_num,])
            Y_weight <- t(as.matrix(output$Y_weight))
            aux <- as.matrix(w_super*Y_weight)
            delta <- as.matrix(output$delta[[j_num]]) %*% aux
            colnames(delta) <- paste0(j, "_", i)
            if(i == 1){
                gene_contribution[[j_num]] <- delta
            }
            if(i > 1){
                gene_contribution[[j_num]] <- cbind(gene_contribution[[j_num]],
                                                    delta)
            }
            j_num <- j_num + 1
        }
    }
    return(list("superpathway_score" = superpathway_score,
                "celltype_contribution" = celltype_contribution,
                "gene_contribution" = gene_contribution))
}
