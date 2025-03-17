#' @title Derive superpathway recapitulation
#'
#' @param model_object A \link{superpathway.fit.model-class} object passed
#' from \link{singISTrecapitulations}
#' @param data_original A matrix with the superpathway's score as returned
#' by \link{derive_contributions} for the non-singIST treated samples,
#' passed from \link{singISTrecapitulations}
#' @param data_singIST A matrix with the superpathway's score as returned
#' by \link{derive_contributions} for the singIST treated samples,
#' passed from \link{singISTrecapitulations}
#'
#' @import checkmate stats
#' @returns
#' An object `data.frame` with the variables: `pathway` name as indicated in
#' `model_object`, `recapitulation` with the superpathway recapitulation
#' @export
#' @examples
#' model <- example_superpathway_fit_model
#' mapped <- example_mapping_organism
#' singIST_samples <- biological_link_function(mapped, model)$singIST_samples
#' original <- derive_contributions(model, singIST_samples)
#' derived <- derive_contributions(model,
#' slot(model, "model_fit")$predictor_block)
#' superpathway_recap(model, original$superpathway_score,
#' derived$superpathway_score)
superpathway_recap <- function(model_object, data_original, data_singIST){
    checkmate::assert_class(model_object, "superpathway.fit.model")
    # Identify indices of base class and target class samples
    base_class <- model_object@superpathway_input@base_class
    target_class <- model_object@superpathway_input@target_class
    indices_base <- which(model_object@superpathway_input@sample_class ==
                            base_class)
    indices_target <- which(model_object@superpathway_input@sample_class ==
                                target_class)
    # Derive reference recapitulation
    Omega <- stats::median(data_original[indices_target]) -
                stats::median(data_original[indices_base])
    # Derive predicted recapitulation
    Omega_prime <- stats::median(data_singIST[indices_base]) -
                    stats::median(data_original[indices_base])
    # Compute predicted recapitulation as a fraction of reference recapitulation
    if(abs(Omega) < .Machine$double.eps){
        message("Reference recapitulation is 0")
        recapitulation <- NULL
    }else{
        recapitulation <- round(100*Omega_prime/Omega, 2)
    }
    superpathway_info <- model_object@superpathway_input@superpathway_info
    output <- data.frame(
        "pathway"= superpathway_info@pathway_info@standard_name,
        "recapitulation" = recapitulation)
    return(output)
}

#' @title Derive cell type recapitulation
#'
#' @param model_object A \link{superpathway.fit.model-class} object passed
#' from \link{singISTrecapitulations}
#' @param data_original A matrix with the cell type contributions as returned
#' by \link{derive_contributions} for the non-singIST treated samples,
#' passed from \link{singISTrecapitulations}
#' @param data_singIST A matrix with the cell type contributions as returned
#' by \link{derive_contributions} for the singIST treated samples,
#' passed from \link{singISTrecapitulations}
#' @import checkmate stats
#'
#' @returns
#' A `data.frame` object with the variables: `pathway` name, `celltype` with
#' the cell type name, `recapitulation` with the cell type recapitulation, and
#' `reference` with the cell type reference recapitulation
#' @export
#' @examples
#' model <- example_superpathway_fit_model
#' mapped <- example_mapping_organism
#' singIST_samples <- biological_link_function(mapped, model)$singIST_samples
#' original <- derive_contributions(model, singIST_samples)
#' derived <- derive_contributions(model,
#' slot(model, "model_fit")$predictor_block)
#' celltype_recap(model, original$celltype_contribution,
#' derived$celltype_contribution)
celltype_recap <- function(model_object, data_original, data_singIST){
    checkmate::assert_class(model_object, "superpathway.fit.model")
    # Identify indices of base class and target class samples
    base_class <- model_object@superpathway_input@base_class
    target_class <- model_object@superpathway_input@target_class
    indices_base <- which(model_object@superpathway_input@sample_class ==
                            base_class)
    indices_target <- which(model_object@superpathway_input@sample_class ==
                                target_class)
    superpathway_info <- model_object@superpathway_input@superpathway_info
    pathway_name <- superpathway_info@pathway_info@standard_name
    # Derive reference and predicted recapitulations
    recapitulation <- data.frame("pathway" = c(), "celltype" = c(),
                                    "recapitulation" = c(),
                                    "reference" = c())
    for(b in seq(1, nrow(data_original))){
        # Reference recapitulation
        Gamma <- stats::median(data_original[b, indices_target])-
                    stats::median(data_original[b, indices_base])
        # Predicted recapitulation
        Gamma_prime <- stats::median(data_singIST[b, indices_base]) -
                        stats::median(data_original[b, indices_base])
        # Predicted recapitulation as a fraction of reference recapitulation
        if(abs(Gamma) < .Machine$double.eps){
            message("Reference recapitulation is 0")
            recapitulation <-
                rbind(recapitulation,
                        data.frame("pathway" = pathway_name,
                        "celltype" = rownames(data_original)[b],
                        "recapitulation" = NULL,
                        "reference" = Gamma))
        }else{
            recapitulation <- rbind(
                recapitulation,data.frame("pathway" = pathway_name,
                "celltype" = rownames(data_original)[b],
                "recapitulation" = 100*Gamma_prime/Gamma,
                "reference" = Gamma))
        }
    }
    return(recapitulation)
}

#' @title Derive gene contribution to cell type recapitulation
#'
#' @param model_object A \link{superpathway.fit.model-class} object passed from
#' \link{singISTrecapitulations}
#' @param data_original A matrix with the gene contributions to superpathway's
#' score as returned by \link{derive_contributions} for the non-singIST treated
#' samples, passed from \link{singISTrecapitulations}
#' @param data_singIST A matrix with the gene contributions to superpathway's
#' score as returned by \link{derive_contributions} for the singIST treated
#' samples, passed from \link{singISTrecapitulations}
#' @param cell_reference A matrix with the cell type recapitulations as
#' returned by \link{celltype_recap}
#' @import checkmate
#' @returns
#' A `data.frame` object with the variables: `pathway` name, `celltype` name,
#' `gene` name, `contribution` gene contribution to cell type recapitulation
#' @export
#' @examples
#' model <- example_superpathway_fit_model
#' mapped <- example_mapping_organism
#' singIST_samples <- biological_link_function(mapped, model)$singIST_samples
#' original <- derive_contributions(model, singIST_samples)
#' derived <- derive_contributions(model,
#' slot(model, "model_fit")$predictor_block)
#' # Derive cell type reference
#' cell <- celltype_recap(model, original$celltype_contribution,
#' derived$celltype_contribution)
#' # Compute gene contributions
#' gene_contrib(model, original$gene_contribution, derived$gene_contribution,
#' cell)
gene_contrib <- function(model_object, data_original,
                        data_singIST, cell_reference){
    checkmate::assert_class(model_object, "superpathway.fit.model")
    # Identify indices of base class
    base_class <- model_object@superpathway_input@base_class
    indices_base <- which(model_object@superpathway_input@sample_class ==
                            base_class)
    superpathway_info <- model_object@superpathway_input@superpathway_info
    pathway_name <- superpathway_info@pathway_info@standard_name
    superpathway_info <- model_object@superpathway_input@superpathway_info
    pathway_name <- superpathway_info@pathway_info@standard_name
    celltypes <- model_object@superpathway_input@superpathway_info@celltypes
    gene_contributions <- data.frame("pathway" = c(), "celltype" = c(),
                                        "gene" = c(), "contribution" = c())
    # Derive gene contributions to cell type recapitulation
    b_num <- 1
    for(b in celltypes){
        delta_prime <- data_singIST[[b_num]]
        delta <- data_original[[b_num]]
        delta_tilde <- delta_prime[,indices_base]-delta[,indices_base]
        # Compute cell type reference recapitulation
        Gamma <- cell_reference[cell_reference$celltype == b, "reference"]
        num_genes <- length(rownames(delta_prime))
        if(abs(Gamma) < .Machine$double.eps){
            contribution <- NULL
        }else{
            contribution <- 100*delta_tilde[,1]/Gamma
        }
        gene_contributions <- rbind(gene_contributions, data.frame(
            "pathway" = rep(pathway_name, num_genes),
            "celltype" = rep(b, num_genes),
            "gene" = rownames(delta_prime),
            "contribution" = contribution
        ))
        b_num <- b_num + 1
    }
    return(gene_contributions)
}

#' @title Compute singIST recapitulations
#'
#' @description
#' This method provides with all singIST recapitulations; superpathway
#' recapitulations; cell type recapitulations; gene contributions to
#' cell type recapitulations. The procedure encompasses the execution of the
#' biological link function, the derivation of the predictor scores
#' (superpathway, cell type and gene scores), and their use to compute
#' the predicted recapitulations as a fraction of the reference recapitulation.
#'
#' @param object A \link{mapping.organism-class} object for which to calculate
#' the recapitulations against the fitted superpathway model
#' @param model_object A \link{superpathway.fit.model-class} object used to
#' calculate the recapitulations
#' @param ... Other parameters to pass onto \link{biological_link_function}
#' @import checkmate
#'
#' @returns
#' A list with; a `data.frame` object with the superpathway recapitulation,
#' containing variables `pathway` name, `recapitulation`, `p_val` with the
#' global significance test of the fitted model as provided in `model_object`,
#' and `target_organism` with the target class of the disease model as provided
#' in `model_object`; a `data.frame` with the cell type recapitulation,
#' containing variables `pathway` name, `celltype`, `recapitulation`,
#' `orthology` with the percentage of observed one-to-one orthology coverage
#' - if all cell types have the same gene set this value is constant -, and
#' `target_organism`; a `data.frame` object with the gene contributions to
#' cell type recapitulation, containing variables `pathway`, `celltype`,
#' `gene` name, `contribution` indicating the gene contribution to cell type
#' recapitulation, and `target_organism`.
#' @export
#' @examples
#' singISTrecapitulations(example_mapping_organism,
#' example_superpathway_fit_model)
singISTrecapitulations <- function(object, model_object, ...){
    checkmate::assert_class(object, "mapping.organism")
    checkmate::assert_class(model_object, "superpathway.fit.model")
    # Derive singIST treated samples
    linkFunction <- biological_link_function(object, model_object, ...)
    C <- model_object@model_fit$predictor_block
    C_prime <- linkFunction$singIST_samples
    # Derive contributions
    contributions_ref <- derive_contributions(model_object, C)
    contributions_singIST <- derive_contributions(model_object, C_prime)
    # Derive superpathway recapitulation
    superpathway <- superpathway_recap(
        model_object, contributions_ref$superpathway_score,
        contributions_singIST$superpathway_score
        )
    # Derive cell type recapitulation
    celltype <- celltype_recap(
        model_object, contributions_ref$celltype_contribution,
        contributions_singIST$celltype_contribution
    )
    # Derive gene contribution to cell type recapitulation
    gene <- gene_contrib(
        model_object, contributions_ref$gene_contribution,
        contributions_singIST$gene_contribution, celltype)
    # Compute observed one-to-one orthology coverage gene set per cell type
    for(b in seq(1, length(object@celltype_mapping))){
        length_set <- length(model_object@model_fit$observed_gene_sets[[b]])
        set <- linkFunction$orthologs[[b]]$output_gene
        observed_set <- sum(set %in% rownames(object@counts))
        coverage <- observed_set/length_set
        # Add observed orthology coverage
        celltype[b, "orthology"] <- 100*coverage
    }
    # Add p value of global significance test
    pval <- model_object@model_validation$pvalue_global_significance
    superpathway[, "p_val"] <- pval
    # Add mapped organism target class to outputs
    superpathway[, "target_organism"] <- object@target_class
    celltype[, "target_organism"] <- object@target_class
    gene[, "target_organism"] <- object@target_class
    celltype[, "reference"] <- NULL
    rownames(gene) <- NULL
    return(list("superpathway" = superpathway, "celltype" = celltype,
                "gene" = gene, "FC" = linkFunction$FC))
}
