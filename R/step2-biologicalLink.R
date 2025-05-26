#' @title Cell type mapping
#' @description
#' For a given \link{mapping.organism-class} it updates the variable
#' `celltype_cluster` so that each element of it is updated accordingly
#' to the mapped cell types as indicated in 
#' `names(slot(object, "celltype_mapping"))`.
#'
#' @param object A \link{mapping.organism-class} object
#'
#' @import checkmate
#' @returns
#' A \link{mapping.organism-class} object with the `slot(object, "counts")`
#' slot updated, for the variable `celltype_cluster` with the cell types
#' according to the mapping defined in `slot(object, "celltype_mapping")`.
#' @export
#' @examples
#' data(example_mapping_organism)
#' data <- example_mapping_organism
#' new_object <- celltype_mapping(data)
#' slot(new_object, "counts")$celltype_cluster
celltype_mapping <- function(object){
    checkmate::assert_class(object, "mapping.organism")
    output <- object
    # Avoid spaces as FindMarkers does not identify them
    names(output@celltype_mapping) <- gsub(" ", "_",
                                            names(output@celltype_mapping))
    # Rename cell types based on cell type mapping
    output@counts$celltype_cluster <- base::unlist(
        base::lapply(output@counts$celltype_cluster, function(x){
            var_name <- names(output@celltype_mapping)[
                vapply(output@celltype_mapping, function(vals) x %in% vals,
                        FUN.VALUE = logical(1))]
            if(length(var_name) > 0) var_name else NA
        }
        )
    )
    # Remove cell types with NA values as those do not have a mapping
    output@counts <- output@counts[, !is.na(output@counts$celltype_cluster)]
    return(output)
}

#' @title Compute differentially expressed genes with FindMarkers
#' @description
#' Computes differentially expressed genes with FindMarkers for the conditions
#' indicated.
#'
#' @param object A \link{mapping.organism-class} object with `Idents(object)`
#' assigned to variables with the conditions being tested.
#' @param condition_1 A vector with the elements of the first factor to perform
#' the hypothesis test. By default the mapped cell types
#' `condition_1 = names(slot(object, "celltype_mapping"))`
#' @param condition_2 A vector with the elements of the second factor to perform
#' the hypothesis test with. By default the class of the organism
#' `condition_2 = c(slot(object, "target_class"), slot(object, "base_class"))`
#' @param logfc.treshold Sets the minimum log-fold change (logFC) cutoff for
#' identifying differentially expressed genes (DEGs). By default
#' `logfc.treshold = 0.25`
#' @param ... Other parameters to pass onto `Seurat::FindMarkers()`
#' @param assay Specific assay being used for analysis. By default
#' `assay = RNA`.
#' @import checkmate Seurat
#' @returns
#' A data frame containing the `Seurat::FindMarkers` output for the interaction
#' of two conditions of the `slot(object, counts)` matrix
#' @export
#' @examples
#' # Set the identities
#' data(example_mapping_organism)
#' data_organism <- example_mapping_organism
#' data <- celltype_mapping(data_organism)
#' slot(data, "counts")$test <- paste0(slot(data, "counts")$celltype_cluster,
#' "_", slot(data, "counts")$class)
#' SeuratObject::Idents(slot(data, "counts")) <- "test"
#' diff_expressed(data)
diff_expressed <- function(object, condition_1 = c(), condition_2 = c(),
                            logfc.treshold = 0.25, assay = "RNA", ...){
    checkmate::assert_class(object, "mapping.organism")
    if(is.null(condition_1)){condition_1 <-
        names(object@celltype_mapping)[lengths(object@celltype_mapping) > 0]}
    if(is.null(condition_2)){condition_2 <-
        c(object@target_class, object@base_class)}
    counts <- object@counts
    # FindMarkers function by row of dataset
    apply_function <- function(row, data = counts) {
        logFC <- Seurat::FindMarkers(
            object = data, ident.1 = row[1], ident.2 = row[2], assay = assay,
            slot = "data", logfc.threshold = logfc.treshold, ...)
        return(logFC)
    }
    # Combinations to test
    combinations <- base::outer(condition_1, condition_2 , paste, sep = "_")
    # Apply the function to each row of the data frame
    output <- base::apply(combinations, 1, apply_function)
    names(output) <- condition_1
    return(output)
}

#' @title Orthology mapping
#' @description
#' Performs the one-to-one orthology mapping between the mapped disease model
#' object \link{mapping.organism-class} to the reference (human) organism
#' of the \link{superpathway.fit.model-class} object.
#'
#' @param object A \link{mapping.organism-class} object
#' @param model_object A \link{superpathway.fit.model-class} object
#' @param from_species A character indicating the reference organism for which
#' the parameter `model_fit` has information from.
#' @param to_species A character indicating the mapped organism for which the
#' parameter `object` has information from. By default `mmusculus`.
#' @param annotation_to_species A character indicating the gene identifier
#' annotation used for the `to_species`. Note this should match with the gene
#' names in `slot(object, counts)`. By default `external_gene_name`. If `NULL`
#' the `annotation_to_species` is inferred with \link{detect_gene_type}, note
#' this might take time.
#' @import biomaRt data.table
#' @returns
#' A list with the gene sets per cell type with the one-to-one orthology
#' @export
#' @examples
#' # Case without stating the gene annotation of the mapping.organisms object
#' # note this will take longer to execute
#' data(example_mapping_organism)
#' data_organism <- example_mapping_organism
#' data(example_superpathway_fit_model)
#' data_model <- example_superpathway_fit_model
#' orthology_mapping(data_organism, data_model, "hsapiens",
#' annotation_to_species = NULL)
#' # Case assuming the gene annotation of the mapping.organism object is
#' # by default "external_gene_name" this is faster
#' orthology_mapping(data_organism, data_model, "hsapiens")
orthology_mapping <- function(object, model_object, from_species,
                                to_species = "mmusculus",
                                annotation_to_species = "external_gene_name"){
    # Connect to Ensembl
    mart_from <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                dataset = paste0(from_species, "_gene_ensembl"))
    mart_to <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                dataset = paste0(to_species, "_gene_ensembl"))
    if(is.null(annotation_to_species)){
        genes_mapped <- rownames(object@counts)
        annotation_to_species <- detect_gene_type(genes_mapped, mart_to)
    }
    gene_sets <- unique(unlist(model_object@model_fit$observed_gene_sets))
    annotation_ref <- detect_gene_type(gene_sets, mart_from)
    # Retrieve ortholog for each observed gene set in the reference species
    gene_set_orthologs <- lapply(
        seq_along(model_object@model_fit$observed_gene_sets),
        function(i, annotation_from = annotation_ref){
            gene_set <- model_object@model_fit$observed_gene_sets[[i]]
            orthologs <- retrieve_one2one_orthologs(
                annotation = annotation_from, gene_set = gene_set,
                mart = mart_from, from_species = from_species,
                to_species = to_species)
            target_gene_symbols <- biomaRt::getBM(
                attributes = c("ensembl_gene_id", annotation_to_species),
                filters = "ensembl_gene_id", values = orthologs$ortholog,
                mart = mart_to)
            data.table::setDT(target_gene_symbols)
            data.table::setnames(target_gene_symbols,
                                new = c("ensembl_gene_id", "output_gene"))
            final_orthologs <- base::merge(
                orthologs, target_gene_symbols, by.x ="ortholog",
                by.y = "ensembl_gene_id", all.x = TRUE)
            return(final_orthologs[, c("output_gene", "input_gene")])
        })
        return(gene_set_orthologs)
}

#' @title Derive singIST treated samples
#'
#' @param object A \link{mapping.organism-class} object passed from
#' \link{biological_link_function}.
#' @param model_object A \link{superpathway.fit.model-class} passed from
#' \link{biological_link_function}
#' @param orthologs A list of `data.table` objects, as returned by
#' \link{orthology_mapping} with the one-to-one orthologs of each gene set per
#' cell type
#' @param logFC A list of `data.frame` objects, as returned by
#' \link{diff_expressed}, with the log Fold Changes computed by FindMarkers
#' for each gene set per cell type
#' @returns
#' A list object with the singIST treated samples predictor block matrix and
#' a list of Fold Changes for each cell type used to compute the singIST
#' treated samples.
#' @export
#' @examples
#' # Orthology mapping
#' data(example_mapping_organism)
#' data_organism <- example_mapping_organism
#' data(example_superpathway_fit_model)
#' data_model <- example_superpathway_fit_model
#' orthologs <- orthology_mapping(data_organism, data_model, "hsapiens")
#' # Set the identities
#' # Cell type mapping
#' data <- celltype_mapping(data_organism)
#' slot(data, "counts")$test <- paste0(slot(data, "counts")$celltype_cluster,
#' "_", slot(data, "counts")$class)
#' SeuratObject::Idents(slot(data, "counts")) <- "test"
#' logFC <- diff_expressed(data)
#' singIST_treat(data_organism, data_model, orthologs, logFC)
singIST_treat <- function(object, model_object, orthologs, logFC){
    samples <- which(model_object@superpathway_input@sample_class ==
                        model_object@superpathway_input@base_class)
    predictor_block <- model_object@model_fit$predictor_block
    cells <- as.vector(which(lengths(object@celltype_mapping) > 0))
    # Update logFC names to remove slashes
    names(logFC) <- gsub("_", " ", names(logFC))
    FC <- vector("list", length(cells))
    for(b in cells){
        genes <- base::intersect(
            orthologs[[b]]$output_gene, rownames(object@counts))
        c <- names(object@celltype_mapping)[b]
        if(length(genes) == 0){
            FC[[c]] <- FC_aux
            next
            }
        FC_aux <- logFC[[c]][rownames(logFC[[c]]) %in% genes, , drop = FALSE]
        significant_genes <- FC_aux[ , "p_val_adj"] <= 0.05
        if(nrow(FC_aux) == 0){
            FC[[c]] <- data.frame()
            next
            }
        if(sum(significant_genes) == 0){
            indices_match <- match(rownames(FC_aux), orthologs[[b]]$output_gene)
            FC_aux[!significant_genes, "avg_log2FC"] <-
                rep(0, sum(!significant_genes))
            rownames(FC_aux) <- paste0(
                c, "*", orthologs[[b]][indices_match, ]$input_gene)
            FC_aux <- FC_aux[, c("avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
            colnames(FC_aux)[1] <- "r_g^b"
            FC[[c]] <- FC_aux
            next
            }
        FC_aux[significant_genes, "avg_log2FC"] <-
            sign(as.numeric(FC_aux[significant_genes, "avg_log2FC"]))*
            2^FC_aux[significant_genes, "avg_log2FC"]
        FC_aux[!significant_genes, "avg_log2FC"] <-
            rep(0, sum(!significant_genes))
        indices_match <- match(rownames(FC_aux), orthologs[[b]]$output_gene)
        rownames(FC_aux) <- paste0(c, "*",
                                    orthologs[[b]][indices_match, ]$input_gene)
        predictor_block <- FCtoExpression(model_object, b, samples,
                                            predictor_block, FC_aux)
        FC_aux <- FC_aux[, c("avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
        colnames(FC_aux)[1] <- "r_g^b"
        FC[[c]] <- FC_aux
    }
    return(list("singIST_samples" = predictor_block, "FC" = FC))
}

#' @title Biological link function
#' @description
#' Maps the organism information in \link{mapping.organism-class} and
#' \link{superpathway.fit.model-class} to obtain the "singIST treated samples"
#' with the simulated human. The biological link function involves the
#' cell type mapping, orthology mapping and fold change computation.
#' @param object A \link{mapping.organism-class} class with the disease model
#' data
#' @param model_object A \link{superpathway.fit.model-class} with the fitted
#' model
#' @param object_gene_identifiers Annotation of gene identifiers used in
#' `object`. By default `external_gene_name`. If `NULL` \link{orthology_mapping}
#' infers the gene identifiers of `object`, note this may add execution time.
#' @param model_species Organism for which `model_object` has been trained. By
#' `default` `hsapiens`.
#' @param ... Other parameters to pass onto \link{diff_expressed} or function
#' \link{orthology_mapping}
#' @import checkmate SeuratObject
#' @returns
#' A list with; ortholog gene sets as returned by \link{orthology_mapping};
#' a list with the Fold Changes used; singIST treated samples as returned by
#' \link{singIST_treat}
#' @export
biological_link_function <- function(
        object, model_object, object_gene_identifiers = "external_gene_name",
        model_species = "hsapiens", ...){
    # Cell type and orthology mapping
    message("Cell type mapping...")
    object <- celltype_mapping(object)
    object@counts$test <- paste0(object@counts$celltype_cluster, "_",
                                    object@counts$class)
    SeuratObject::Idents(object@counts) <- "test"
    message("Computing Fold Changes with FindMarkers...")
    logFC <- diff_expressed(object, ...)
    message("Orthology mapping...")
    to_species <- paste0(tolower(substr(object@organism,1,1)),
                            tolower(sub(".* ", "", object@organism)))
    # Remove "_" from cell type name once diff_expressed is executed
    names(object@celltype_mapping) <- gsub("_", " ",
                                            names(object@celltype_mapping))
    if(to_species != model_species){
        orthologs <- orthology_mapping(
            object, model_object, to_species = to_species,
            annotation_to_species = object_gene_identifiers, 
            from_species = model_species, ...)
    }else{ # Case where no orthology mapping should be applied 
        orthologs <-lapply(seq_along(model_object@model_fit$observed_gene_sets),
                        function(i){
                        sets <- model_object@model_fit$observed_gene_sets[[i]]
                        aux <- data.table("input_gene" = sets,
                                            "output_gene" = sets)
                        aux
                        })
    }
    # singIST treated samples
    message("Deriving singIST treated samples...")
    singIST_samples <- singIST_treat(object, model_object, orthologs, logFC)
    # Set names
    names(orthologs) <- names(object@celltype_mapping)
    return(list("orthologs" = orthologs,
            "singIST_samples" = singIST_samples$singIST_samples,
            "FC" = singIST_samples$FC))
}
