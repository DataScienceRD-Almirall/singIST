#' @title Check if parameter format is consistent
#'
#' @description
#' For the wrapper functions \link{multiple_fitOptimal} and
#' \link{multiple_singISTrecapitulations} one must pass multiple parameters. To
#' check for the consistency of such parameters we use this function with the
#' logic; if the parameter passed is `NULL` or its length is 1, it is assumed
#' that the desired list of parameters is the repetition of such; if the
#' parameter passed is a vector whose length is the number of objects that
#' the wrappers iterate on, then the function returns a list whose elements
#' are each of the vector elements; otherwise if the parameters are a vector
#' whose length does not match with the number of objects to iterate on
#' then the function stops
#'
#' @param parameter The parameters passed from either
#' \link{multiple_fitOptimal} or \link{multiple_singISTrecapitulations}
#' @param objectLength The number of objects that the wrapper functions
#' iterate on
#'
#' @returns
#' A list with the repetition of the parameter
#'
#' @import checkmate
#' @export
#' @examples
#' # NULL parameter case
#' parameter <- NULL
#' objectLength <- 10
#' multiple_check(parameter, objectLength)
#'
#' # Parameter equal for all elements
#' parameter <- FALSE
#' objectLength <- 5
#' multiple_check(parameter, objectLength)
#'
#' # Parameter differing for all elements
#' parameter <- c(1, 6, 7, 8, 9)
#' objectLength <- 5
#' multiple_check(parameter, objectLength)
multiple_check <- function(parameter, objectLength){
    output <- vector("list", length = objectLength)
    # Check if parameter is null then all the list should be NULL
    if(is.null(parameter)){
        for(i in seq(1, objectLength)){
            output <- vector("list", length = objectLength)
        }
    }else{
        if(length(parameter) == 1){
            for(i in seq(1, objectLength)){
                output[[i]] <- parameter
            }
        }else if(length(parameter) == objectLength){
            for(i in seq(1, objectLength)){
                output[[i]] <- parameter[[i]]
            }
        }else{
            stop("Parameters should either have length 1 (equal for all objects
                to assess) or same length as number of objects.")
        }
    }
    return(output)
}

#' @title Multiple Cross validation and fit of asmbPLSDA
#'
#' @description
#' Use \link{fitOptimal} for multiple \link{superpathway.input-class}
#' objects. This wrapper is useful if one wants to assess multiple superpathways
#' for analyses and needs to train its respective optimal models.
#'
#' @param object A list whose elements are \link{superpathway.input-class}
#' objects to use \link{fitOptimal} to
#' @param parallel A vector whose elements are `parallel` parameters for each
#' object as requested by \link{fitOptimal}
#' @param measure A vector whose elements are `measure` parameters for each
#' object as requested by \link{fitOptimal}
#' @param expected_measure_increase A vector whose elements are
#' `expected_measure_increase` parameters for each object as requested by
#' \link{fitOptimal}
#' @param maxiter A vector whose elements are `maxiter` parameters for each
#' object as requested by \link{fitOptimal}
#' @param global_significance_full A vector whose elements are
#' `global_significance_full` parameters for each object as requested by
#' \link{fitOptimal}
#' @param CIP.GIP_significance_full A vector whose elements are
#' `CIP.GIP_significance_full` parameters for each object as requested by
#' \link{fitOptimal}
#' @param npermut A vector whose elements are `npermut` parameters for each
#' object as requested by \link{fitOptimal}
#' @param nbObsPermut A vector whose elements are `nbObsPermut` parameters for
#' each object as requested by \link{fitOptimal}
#' @param type A vector whose elements are `type` parameters for
#' each object as requested by \link{fitOptimal}
#' @param nsubsampling A vector whose elements are `nsubsampling` parameters for
#' each object as requested by \link{fitOptimal}
#' @param Method A vector whose elements are `Method` parameters for
#' each object as requested by \link{fitOptimal}
#'
#' @import checkmate
#' @returns
#' A list of \link{superpathway.fit.model-class} objects
#' @export
#' @examples
#' data(example_superpathway_input)
#' data <- example_superpathway_input
#' models <- list(data, data)
#' # Example with different options
#' multiple_model <- multiple_fitOptimal(models, type = c("jackknife",
#' "subsampling"), nsubsampling = c(NULL, 10), npermut = c(10,15))
multiple_fitOptimal <- function(
        object = list(), parallel = c(FALSE),
        measure = c("B_accuracy"), expected_measure_increase = c(0.005),
        maxiter = c(100), global_significance_full= c(FALSE),
        CIP.GIP_significance_full = c(FALSE), npermut = c(100),
        nbObsPermut = c(NULL), type = c("jackknife"),
        nsubsampling = c(100), Method = c(NULL)
        ){
    # Check that multiple objects are provided in the proper format
    nobjects <- length(object)
    checkmate::assert_list(object, any.missing = FALSE, all.missing = FALSE)
    checkmate::assert_true(nobjects > 1)
    for(element in object){
        checkmate::assert_class(element, "superpathway.input")
    }
    # If parameter provided is NULL then the rest is NULL otherwise provide all
    parallel <- multiple_check(parallel, nobjects)
    measure <- multiple_check(measure, nobjects)
    expected_measure_increase <- multiple_check(
        expected_measure_increase, nobjects)
    maxiter <- multiple_check(maxiter, nobjects)
    global_significance_full<- multiple_check(
        global_significance_full, nobjects)
    CIP.GIP_significance_full <- multiple_check(
        CIP.GIP_significance_full, nobjects)
    npermut <- multiple_check(npermut, nobjects)
    nbObsPermut <- multiple_check(nbObsPermut, nobjects)
    type <- multiple_check(type, nobjects)
    nsubsampling <- multiple_check(nsubsampling, nobjects)
    # fitOptimal for each object in the list
    model <- vector("list", length = nobjects)
    for(l in seq(1, nobjects)){
        pathway <- object[[l]]@superpathway_info@pathway_info@standard_name
        names(model)[[l]] <- pathway
        message("Object ", pathway)
        model[[l]] <- fitOptimal(
            object = object[[l]], parallel = parallel[[l]],
            measure = measure[[l]], expected_measure_increase =
            expected_measure_increase[[l]], maxiter = maxiter[[l]],
            global_significance_full= global_significance_full[[l]],
            CIP.GIP_significance_full = CIP.GIP_significance_full[[l]],
            npermut = npermut[[l]], nbObsPermut = nbObsPermut[[l]],
            type = type[[l]], nsubsampling = nsubsampling[[l]],
            Method = Method[[l]]
            )
    }
    return(model)
}

#' @title Compute singIST recapitulations for multiple superpathways
#'
#' @description
#' Use \link{singISTrecapitulations} for multiple
#' \link{superpathway.input-class} objects against the same
#' \link{mapping.organism-class} object. This wrapper is useful if one wants
#' to assess multiple superpathways against the same
#' \link{mapping.organism-class}
#'
#' @param object A \link{mapping.organism-class} object
#' @param model_object A list whose elements
#' are \link{superpathway.fit.model-class}
#' @param from_species A list of characters indicating the organism of each
#' `model_object` element. By default `list("hsapiens")` which assumes the same
#' organism across all elements of `model_object` parameter
#' @param ...  Other parameters to pass onto \link{biological_link_function}
#'
#' @import checkmate
#' @returns
#' A list with the row binded `data.frame` for each superpathway assessed for
#' the superpathway and cell type recapitulations, and gene contributions
#' to the former.
#' @export
#' @examples
#' data(example_superpathway_input)
#' data_model <- example_superpathway_input
#' models <- list(data_model, data_model)
#' # Example with different options
#' multiple_model <- multiple_fitOptimal(models, type = c("jackknife",
#' "subsampling"), nsubsampling = c(NULL, 10), npermut = c(10,15))
#' data(example_mapping_organism)
#' data_organism <- example_mapping_organism
#' multiple_singISTrecapitulations(data_organism, multiple_model,
#' from_species = list("hsapiens", "hsapiens"))
multiple_singISTrecapitulations <- function(
        object, model_object = list(), from_species = list("hsapiens"), ...){
    # Check that multiple objects are provided in the proper format
    nmodels <- length(model_object)
    checkmate::assert_list(model_object, any.missing = FALSE,
                            all.missing = FALSE)
    checkmate::assert_true(nmodels > 1)
    for(element in model_object){
        checkmate::assert_class(element, "superpathway.fit.model")
    }
    from_species <- multiple_check(from_species, nmodels)
    # Initialize output
    aux_superpathway <- vector("list", length = nmodels)
    aux_celltype <- vector("list", length = nmodels)
    aux_gene <- vector("list", length = nmodels)
    aux_FC <- vector("list", length = nmodels)
    aux_orthologs <- vector("list", length = nmodels)
    for(i in seq(1, nmodels)){
        pathway <- model_object[[i]]@superpathway_input@superpathway_info
        message("Executing superpathway ", pathway@pathway_info@standard_name)
        aux <- singISTrecapitulations(
            object, model_object[[i]], from_species = from_species[[i]], ...)
        aux_superpathway[[i]] <- aux$superpathway
        aux_celltype[[i]] <- aux$celltype
        aux_gene[[i]] <- aux$gene
        aux_FC[[i]] <- aux$FC
        names(aux_FC)[i] <- pathway@pathway_info@standard_name
        aux_orthologs[[i]] <- aux$orthologs
    }
    output <- vector("list", length = 5)
    names(output) <- c("superpathway", "celltype", "gene", "FC", "orthologs")
    # Aggregate output in a single data.frame
    output$superpathway <- base::do.call(rbind, aux_superpathway)
    output$celltype <- base::do.call(rbind, aux_celltype)
    output$gene <- base::do.call(rbind, aux_gene)
    output$orthologs <- aux_orthologs
    for(i in seq(1, nmodels)){
        output$FC[[i]] <- aux_FC[[i]]
        names(output$FC)[i] <- names(aux_FC)[i]
    }
    return(output)
}

#' @title Render multiple singISTrecapitulation outputs
#'
#' @description
#' Render output of \link{multiple_singISTrecapitulations} for multiple
#' disease models and superpathways, the output is friendly for visualizing
#' the results
#' @param objects A list as retuned by \link{multiple_singISTrecapitulations}
#' @import checkmate purrr
#' @returns
#' A list with the row binded `data.frame` for each superpathway assessed for
#' the superpathway and cell type recapitulations, gene contributions
#' to the former, and fold changes. These row binds are performed for all
#' disease models and superpathways.
#' @export
render_multiple_outputs <- function(objects = list()){
    checkmate::assert_true(length(objects) >= 2)
    superpathways <- do.call(rbind, lapply(seq_along(objects), function(i){
        objects[[i]]$superpathway
        })
        )
    celltypes <- do.call(rbind, lapply(seq_along(objects), function(i){
        objects[[i]]$celltype
        })
        )
    genes <- do.call(rbind, lapply(seq_along(objects), function(i){
        orthologs <- purrr::pmap_int(
            objects[[i]]$gene[, c("gene", "pathway", "celltype")],
            function(gene, pathway, celltype){
                names(objects[[i]]$orthologs) <- names(objects[[i]]$FC)
                genes_in_pathway <- 
                    objects[[i]]$orthologs[[pathway]][[celltype]]$input_gene
                if(gene %in% genes_in_pathway) return(1)
                return(0)
            })
        objects[[i]]$gene$orthology <- orthologs
        objects[[i]]$gene
        })
        )
    fc <- do.call(rbind, lapply(seq_along(objects), function(i){# For each model
        all <- do.call(
            rbind,
            lapply(seq_along(objects[[i]]$FC),function(j,data =objects[[i]]$FC){
                combined <- do.call(rbind, data[[j]])
                clean_name <- sub("^[^.]+\\.", "", rownames(combined))
                celltype_name <- sub("\\*.*$", "", clean_name)
                gene_name <- sub("^.*\\*", "", clean_name)
                pathway_name_col <- data.frame("pathway" = rep(names(data)[j],
                                                                nrow(combined)))
                result <- cbind(pathway_name_col, "celltype" = celltype_name,
                                "gene" = gene_name, combined)
                rownames(result) <- NULL
                return(result)
            }))
        all <- cbind(all, "target_organism" =
                        rep(unique(objects[[i]]$superpathway$target_organism),
                            nrow(all)))
        all
    }))
    output <- list("superpathway" = superpathways, "celltype" = celltypes,
                    "gene" = genes, "FC" = fc)
    return(output)
}