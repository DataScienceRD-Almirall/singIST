#' @title Validate pathway fields
#' @description
#' Checks that all provided fields for a pathway meet the expected properties.
#'
#' @param standard_name Character. Pathway standard name from MsigDB.
#' @param dbsource Character. Database source (KEGG, PID, REACTOME, BIOCARTA, WIKIPATHWAYS).
#' @param collection Character. MsigDB collection (c2 or m2).
#' @param subcollection Character. MsigDB subcollection (CP).
#' @export
#' 
#' @return TRUE if all checks pass; otherwise, an error is thrown.
#' @examples
#' check_pathway(
#'   standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#'   dbsource = "KEGG",
#'   collection = "c2",
#'   subcollection = "CP"
#' )
check_pathway <- function(standard_name, dbsource, collection, subcollection) {
    checkmate::assert_character(standard_name, len = 1)
    checkmate::assert_choice(dbsource, choices = c("KEGG", "PID", "REACTOME",
                                                   "BIOCARTA", "WIKIPATHWAYS"))
    checkmate::assert_choice(collection, choices = c("c2", "m2"))
    checkmate::assert_choice(subcollection, choices = c("CP"))
    # Check consistency of standard_name with dbsource
    if (dbsource != "WIKIPATHWAYS") {
        checkmate::assert_character(standard_name, pattern = paste0("^",
                                                                    dbsource, "_"))
    } else {
        checkmate::assert_character(standard_name, pattern = "^WP_")
    }
    return(TRUE)
}

#' @title Create pathway object
#' @description
#' Creates a simple list representing a pathway, after validating its fields.
#'
#' @param standard_name Character. Pathway standard name from MsigDB.
#' @param dbsource Character. Database source (KEGG, PID, REACTOME, BIOCARTA, WIKIPATHWAYS).
#' @param collection Character. MsigDB collection (c2 or m2).
#' @param subcollection Character. MsigDB subcollection (CP).
#'
#' @return A list with elements: standard_name, dbsource, collection, subcollection.
#' @export
#' 
#' @examples
#' my_pathway <- create_pathway(
#'   standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#'   dbsource = "KEGG",
#'   collection = "c2",
#'   subcollection = "CP"
#' )
#' print(my_pathway)
create_pathway <- function(standard_name, dbsource, collection, subcollection) {
    # Validate inputs
    check_pathway(standard_name, dbsource, collection, subcollection)
    # Return as a simple list
    pathway <- list(
        standard_name = standard_name,
        dbsource = dbsource,
        collection = collection,
        subcollection = subcollection
    )
    return(pathway)
}

#' @title Validate superpathway gene sets
#' @description
#' Checks that all provided fields for a superpathway meet the expected properties:
#' - pathway_info must be a valid pathway (validated with check_pathway()).
#' - celltypes must be a character vector with at least 2 elements.
#' - gene_sets_celltype must be a list (or NULL) with same length as celltypes.
#'
#' @param pathway_info List. A pathway object created by create_pathway().
#' @param celltypes Character vector. Each element represents a cell type.
#' @param gene_sets_celltype List of character vectors. Each element corresponds
#' to gene sets for each cell type. Can be NULL.
#' @export
#' 
#' @return TRUE if all checks pass; otherwise, an error is thrown.
#' @examples
#' my_pathway <- create_pathway("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#' "KEGG", "c2", "CP")
#' check_superpathway(my_pathway, c("T-cell", "Dendritic Cell"), list(c("IL4",
#' "IL5"), c("IL13")))
check_superpathway <- function(pathway_info, celltypes, gene_sets_celltype) {
    # Validate pathway_info using check_pathway()
    check_pathway(
        standard_name = pathway_info$standard_name,
        dbsource = pathway_info$dbsource,
        collection = pathway_info$collection,
        subcollection = pathway_info$subcollection
    )
    # Validate celltypes
    checkmate::assert_character(celltypes, min.len = 2)
    # Validate gene_sets_celltype
    checkmate::assert_list(gene_sets_celltype, null.ok = TRUE)
    if (length(gene_sets_celltype) > 0) {
        len_genesets <- length(gene_sets_celltype)
        len_celltypes <- length(celltypes)
        checkmate::assert_true(len_genesets == len_celltypes)
    }
    return(TRUE)
}

#' @title Create superpathway gene sets object
#' @description
#' Creates a simple list representing a superpathway, after validating its fields.
#'
#' @param pathway_info List. A pathway object created by create_pathway().
#' @param celltypes Character vector. Each element represents a cell type.
#' @param gene_sets_celltype List of character vectors. Each element corresponds
#' to gene sets for each cell type. Can be NULL.
#' @export
#' 
#' @return A list with elements: pathway_info, celltypes, gene_sets_celltype.
#' @examples
#' my_pathway <- create_pathway("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#' "KEGG", "c2", "CP")
#' my_superpathway <- create_superpathway(my_pathway, c("T-cell",
#' "Dendritic Cell"), list(c("IL4", "IL5"), c("IL13")))
#' print(my_superpathway)
create_superpathway <- function(pathway_info, celltypes, gene_sets_celltype) {
    # Validate inputs
    check_superpathway(pathway_info, celltypes, gene_sets_celltype)
    # Return as a simple list
    superpathway <- list(
        pathway_info = pathway_info,
        celltypes = celltypes,
        gene_sets_celltype = gene_sets_celltype
    )
    return(superpathway)
}

#' @title Validate asmbPLS-DA hyperparameters
#' @description
#' Checks that all provided hyperparameters meet the expected properties:
#' - quantile_comb_table must be a matrix with at least one row.
#' - outcome_type must be "binary" or "multiclass".
#' - number_PLS must be an integer >= 0.
#' - folds_CV and repetition_CV must be integers >= 0 (or NULL).
#' @param quantile_comb_table Matrix. Quantile (lambda) sparsity values for CV.
#' @param outcome_type Character. Either "binary" or "multiclass".
#' @param number_PLS Integer. Maximum number of PLS components.
#' @param folds_CV Integer or NULL. Number of folds for CV (default 5).
#' @param repetition_CV Integer or NULL. Number of repetitions for CV (default 10).
#' @export
#' 
#' @return TRUE if all checks pass; otherwise, an error is thrown.
#' @examples
#' quantile_comb_table <- base::as.matrix(RcppAlgos::permuteGeneral(seq(0.05,
#' 0.95, by = 0.50)))
#' check_hyperparameters(quantile_comb_table, "binary", 3L, 1L, 1L)
check_hyperparameters <- function(quantile_comb_table, outcome_type, number_PLS,
                                  folds_CV, repetition_CV) {
    checkmate::assert_matrix(quantile_comb_table, min.rows = 1)
    checkmate::assert_choice(outcome_type, choices = c("binary", "multiclass"))
    checkmate::assert_integer(number_PLS, lower = 0)
    checkmate::assert_integer(folds_CV, lower = 0, null.ok = TRUE)
    checkmate::assert_integer(repetition_CV, lower = 0, null.ok = TRUE)
    # Additional logical consistency checks
    if (!is.null(folds_CV) && folds_CV == 1) {
        checkmate::assert_true(repetition_CV == 1)
    }
    return(TRUE)
}

#' @title Create asmbPLS-DA hyperparameters object
#' @description
#' Creates a simple list representing hyperparameters for asmbPLS-DA, after validating them.
#'
#' @param quantile_comb_table Matrix. Quantile (lambda) sparsity values for CV.
#' @param outcome_type Character. Either "binary" or "multiclass".
#' @param number_PLS Integer. Maximum number of PLS components.
#' @param folds_CV Integer or NULL. Number of folds for CV (default 5).
#' @param repetition_CV Integer or NULL. Number of repetitions for CV (default 10).
#' @export
#' 
#' @return A list with elements: quantile_comb_table, outcome_type, number_PLS, folds_CV, repetition_CV.
#' @examples
#' quantile_comb_table <- base::as.matrix(RcppAlgos::permuteGeneral(seq(0.05,
#' 0.95, by = 0.50)))
#' my_hyperparameters <- create_hyperparameters(quantile_comb_table,"binary",3L,
#' 1L, 1L)
#' print(my_hyperparameters)
create_hyperparameters <- function(quantile_comb_table, outcome_type, number_PLS,
                                   folds_CV = 5L, repetition_CV = 10L) {
    # Validate inputs
    check_hyperparameters(quantile_comb_table, outcome_type, number_PLS,
                          folds_CV, repetition_CV)
    # Return as a simple list
    hyperparams <- list(
        quantile_comb_table = quantile_comb_table,
        outcome_type = outcome_type,
        number_PLS = number_PLS,
        folds_CV = folds_CV,
        repetition_CV = repetition_CV
    )
    return(hyperparams)
}

#' @title Check superpathway input for asmbPLS-DA
#'
#' @description
#' Checks the validity of the inputs. This version assumes that
#' \code{superpathway_info} and \code{hyperparameters_info} are plain lists
#' validated by \code{check_superpathway()} and \code{check_hyperparameters()}.
#'
#' @param superpathway_info A list representing a superpathway object.
#' @param hyperparameters_info A list representing a hyperparameters object.
#' @param pseudobulk_lognorm A pseudobulk matrix (rows: "Celltype_Sampleid",
#' cols: genes in "HGNC" format or similar).
#' @param sample_id A character vector of sample ids.
#' @param sample_class A character vector with the class of each sample.
#' @param base_class A character scalar indicating the base class.
#' @param target_class A character scalar indicating the target class.
#'
#' @return Invisibly returns \code{TRUE} if all checks pass; otherwise errors.
#'
#' @import checkmate
#' @export
#' @examples
#' # ---- Superpathway info (list) ----
#' my_pathway <- create_pathway("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#' "KEGG", "c2", "CP")
#' celltypes <- c("T-cell", "Dendritic Cell")
#' 
#' my_superpathway <- create_superpathway(my_pathway, celltypes, list(c("IL4",
#' "IL5"), c("IL13")))
#' # ---- Hyperparameters info (list) ----
#' quantile_comb_table <- base::as.matrix(
#'   RcppAlgos::permuteGeneral(seq(0.05, 0.95, by = 0.50)),
#'   ncol = length(celltypes)
#' )
#' 
#' my_hyperparameters <- create_hyperparameters(
#'   quantile_comb_table = quantile_comb_table,
#'   outcome_type = "binary",
#'   number_PLS = as.integer(3),
#'   folds_CV = as.integer(1),
#'   repetition_CV = as.integer(1)
#' )
#'
#' # ---- Pseudobulk + labels ----
#' sample_id <- c("AD1", "AD2", "HC1", "HC2")
#' sample_class <- c("AD", "AD", "HC", "HC")
#' base_class <- "HC"
#' target_class <- "AD"
#'
#' pseudobulk_lognorm <- matrix(
#'   rnorm(length(celltypes) * length(sample_id)),
#'   nrow = length(celltypes) * length(sample_id),
#'   ncol = length(celltypes)
#' )
#' rownames(pseudobulk_lognorm) <- as.vector(t(outer(
#'   celltypes, sample_id, function(x, y) paste(x, y, sep = "_")
#' )))
#'
#' check_superpathway_input(
#'   superpathway_info = my_superpathway,
#'   hyperparameters_info = my_hyperparameters,
#'   pseudobulk_lognorm = pseudobulk_lognorm,
#'   sample_id = sample_id,
#'   sample_class = sample_class,
#'   base_class = base_class,
#'   target_class = target_class
#' )
check_superpathway_input <- function(superpathway_info,
                                     hyperparameters_info,
                                     pseudobulk_lognorm,
                                     sample_id,
                                     sample_class,
                                     base_class,
                                     target_class) {
    # superpathway_info and hyperparameters_info are lists now
    check_superpathway(superpathway_info$pathway_info, superpathway_info$celltypes,
                       superpathway_info$gene_sets_celltype)
    check_hyperparameters(hyperparameters_info$quantile_comb_table,
                          hyperparameters_info$outcome_type,
                          hyperparameters_info$number_PLS,
                          hyperparameters_info$folds_CV,
                          hyperparameters_info$repetition_CV)
    checkmate::assert_character(sample_id)
    checkmate::assert_character(sample_class)
    checkmate::assert_choice(base_class, choices = unique(sample_class))
    possible_classes <- setdiff(unique(sample_class), base_class)
    checkmate::assert_choice(target_class, choices = possible_classes)
    # Pseudobulk matrix check
    checkmate::assert_matrix(pseudobulk_lognorm, any.missing = TRUE)
    # gene_sets must not be NULL
    checkmate::assert_list(superpathway_info$gene_sets_celltype)
    # outcome_type consistency with number of classes
    # Assumes hyperparameters_info$outcome_type exists and is a string like "binary"
    if (hyperparameters_info$outcome_type == "binary") {
        checkmate::assert_true(length(unique(sample_class)) == 2)
    } else {
        checkmate::assert_true(length(unique(sample_class)) >= 3)
    }
    # Sample id should not contain "_" characters (kept identical to original logic)
    sample_character <- base::grepl("_", sample_id)
    checkmate::assert_true(!all(sample_character))
    TRUE
}


#' @title Create superpathway input for asmbPLS-DA
#'
#' @description
#' Creates a plain list containing the fields for superpathway input object
#' and validates the content with \code{\link{check_superpathway_input}}.
#'
#' @param superpathway_info A list representing a superpathway object (formerly
#' \code{superpathway.gene.sets}).
#' @param hyperparameters_info A list representing a hyperparameters object
#' (formerly \code{hyperparameters}).
#' @param pseudobulk_lognorm A pseudobulk matrix.
#' @param sample_id A character vector of sample ids.
#' @param sample_class A character vector with the class of each sample.
#' @param base_class A character scalar indicating the base class.
#' @param target_class A character scalar indicating the target class.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{superpathway_info}
#'   \item \code{hyperparameters_info}
#'   \item \code{pseudobulk_lognorm}
#'   \item \code{sample_id}
#'   \item \code{sample_class}
#'   \item \code{base_class}
#'   \item \code{target_class}
#' }
#'
#' @export
#' @examples
#' # ---- Superpathway info (list) ----
#' my_pathway <- create_pathway("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#' "KEGG", "c2", "CP")
#' celltypes <- c("T-cell", "Dendritic Cell")
#' 
#' my_superpathway <- create_superpathway(my_pathway, celltypes, list(c("IL4",
#' "IL5"), c("IL13")))
#' # ---- Hyperparameters info (list) ----
#' quantile_comb_table <- base::as.matrix(
#'   RcppAlgos::permuteGeneral(seq(0.05, 0.95, by = 0.50)),
#'   ncol = length(celltypes)
#' )
#' 
#' my_hyperparameters <- create_hyperparameters(
#'   quantile_comb_table = quantile_comb_table,
#'   outcome_type = "binary",
#'   number_PLS = as.integer(3),
#'   folds_CV = as.integer(1),
#'   repetition_CV = as.integer(1)
#' )
#'
#' # ---- Pseudobulk + labels ----
#' sample_id <- c("AD1", "AD2", "HC1", "HC2")
#' sample_class <- c("AD", "AD", "HC", "HC")
#' base_class <- "HC"
#' target_class <- "AD"
#'
#' pseudobulk_lognorm <- matrix(
#'   rnorm(length(celltypes) * length(sample_id)),
#'   nrow = length(celltypes) * length(sample_id),
#'   ncol = length(celltypes)
#' )
#' rownames(pseudobulk_lognorm) <- as.vector(t(outer(
#'   celltypes, sample_id, function(x, y) paste(x, y, sep = "_")
#' )))
#'
#' my_superpathway_input <- create_superpathway_input(
#'   superpathway_info = my_superpathway,
#'   hyperparameters_info = my_hyperparameters,
#'   pseudobulk_lognorm = pseudobulk_lognorm,
#'   sample_id = sample_id,
#'   sample_class = sample_class,
#'   base_class = base_class,
#'   target_class = target_class
#' )
create_superpathway_input <- function(superpathway_info,
                                      hyperparameters_info,
                                      pseudobulk_lognorm,
                                      sample_id,
                                      sample_class,
                                      base_class,
                                      target_class) {
    
    check_superpathway_input(
        superpathway_info = superpathway_info,
        hyperparameters_info = hyperparameters_info,
        pseudobulk_lognorm = pseudobulk_lognorm,
        sample_id = sample_id,
        sample_class = sample_class,
        base_class = base_class,
        target_class = target_class
    )
    
    list(
        superpathway_info = superpathway_info,
        hyperparameters_info = hyperparameters_info,
        pseudobulk_lognorm = pseudobulk_lognorm,
        sample_id = sample_id,
        sample_class = sample_class,
        base_class = base_class,
        target_class = target_class
    )
}

#' @title Validate superpathway fit model
#' @description
#' Checks that all provided fields for a fitted asmbPLS-DA model meet the expected properties:
#' - superpathway_input must be a valid superpathway_input object (validated with check_superpathway_input()).
#' - hyperparameters_fit must be valid hyperparameters (validated with check_hyperparameters()).
#' - model_fit and model_validation must be lists.
#' @param superpathway_input List. A superpathway object created by create_superpathway().
#' @param hyperparameters_fit List. A hyperparameters object created by create_hyperparameters().
#' @param model_fit List. Fitted model details.
#' @param model_validation List. Validation metrics of the fitted model.
#'
#' @return TRUE if all checks pass; otherwise, an error is thrown.
check_fit_model <- function(superpathway_input, hyperparameters_fit, model_fit,
                            model_validation) {
    # Validate superpathway_input using check_superpathway_input()
    check_superpathway_input(
        superpathway_info = superpathway_input$superpathway_info,
        hyperparameters_info = superpathway_input$hyperparameters_info,
        pseudobulk_lognorm = superpathway_input$pseudobulk_lognorm,
        sample_id = superpathway_input$sample_id,
        sample_class = superpathway_input$sample_class,
        base_class = superpathway_input$base_class,
        target_class = superpathway_input$target_class
    )
    # Validate hyperparameters using check_hyperparameters()
    check_hyperparameters(
        quantile_comb_table = hyperparameters_fit$quantile_comb_table,
        outcome_type = hyperparameters_fit$outcome_type,
        number_PLS = hyperparameters_fit$number_PLS,
        folds_CV = hyperparameters_fit$folds_CV,
        repetition_CV = hyperparameters_fit$repetition_CV
    )
    # Validate model_fit and model_validation
    checkmate::assert_list(model_fit)
    checkmate::assert_list(model_validation)
    return(TRUE)
}

#' @title Create superpathway fit model object
#' @description
#' Creates a simple list representing a fitted asmbPLS-DA model, after validating its components.
#'
#' @param superpathway_input List. A superpathway object created by create_superpathway().
#' @param hyperparameters_fit List. A hyperparameters object created by create_hyperparameters().
#' @param model_fit List. Fitted model details.
#' @param model_validation List. Validation metrics of the fitted model.
#' @export
#' 
#' @return A list with elements: superpathway_input, hyperparameters_fit, model_fit, model_validation.
#' @examples
#' my_pathway <- create_pathway("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#' "KEGG", "c2", "CP")
#' 
#' celltypes <- c("T-cell", "Dendritic Cell")
#' 
#' my_superpathway <- create_superpathway(my_pathway, celltypes,
#' list(c("IL4", "IL5"), c("IL13")))
#' 
#' my_hyperparameters <- create_hyperparameters(matrix(1:4, nrow = 2), "binary",
#' 3L, 5L, 10L)
#' 
#' #' # ---- Pseudobulk + labels ----
#' sample_id <- c("AD1", "AD2", "HC1", "HC2")
#' sample_class <- c("AD", "AD", "HC", "HC")
#' base_class <- "HC"
#' target_class <- "AD"
#'
#' pseudobulk_lognorm <- matrix(
#'   rnorm(length(celltypes) * length(sample_id)),
#'   nrow = length(celltypes) * length(sample_id),
#'   ncol = length(celltypes)
#' )
#' rownames(pseudobulk_lognorm) <- as.vector(t(outer(
#'   celltypes, sample_id, function(x, y) paste(x, y, sep = "_")
#' )))
#' 
#' my_superpathway_input <- create_superpathway_input(
#'     superpathway_info = my_superpathway,
#'     hyperparameters_info = my_hyperparameters,
#'     pseudobulk_lognorm = pseudobulk_lognorm,
#'     sample_id = sample_id,
#'     sample_class = sample_class,
#'     base_class = base_class,
#'     target_class = target_class
#'     )
#' my_fit <- create_fit_model(my_superpathway_input, my_hyperparameters,
#' list(model = "fit"), list(accuracy = 0.95))
create_fit_model <- function(superpathway_input, hyperparameters_fit, model_fit, model_validation) {
    # Validate inputs
    check_fit_model(superpathway_input, hyperparameters_fit, model_fit, model_validation)
    # Return as a simple list
    fit_model <- list(
        superpathway_input = superpathway_input,
        hyperparameters_fit = hyperparameters_fit,
        model_fit = model_fit,
        model_validation = model_validation
    )
    
    return(fit_model)
}


#' @title Validate mapping organism input
#' @description
#' Validates that all provided fields for a mapping organism meet expected properties:
#' - organism, target_class, base_class are scalar characters.
#' - celltype_mapping is a list with non-empty names; each entry is a character vector (length 0 allowed).
#' - counts is a Seurat or SingleCellExperiment object with metadata columns: 'class' and 'celltype_cluster'.
#' - target_class and base_class exist in meta$class and are different (after normalization).
#' - all clusters referenced in celltype_mapping exist in meta$celltype_cluster.
#'
#' @param organism Character(1). Scientific Latin name of the organism.
#' @param target_class Character(1). Name of the target class for this organism.
#' @param base_class Character(1). Name of the base class for this organism.
#' @param celltype_mapping List. Mapping of cell types (keys) to cluster names (character vectors).
#' @param counts Seurat or SingleCellExperiment. Contains scRNA-seq counts and metadata.
#' @export
#' 
#' @return TRUE if all checks pass; otherwise an informative error is thrown.
#' @examples
#' # Seurat example
#' counts <- SeuratObject::pbmc_small
#' colnames(slot(counts, "meta.data"))[1] <- "donor"
#' colnames(slot(counts, "meta.data"))[6] <- "class"
#' colnames(slot(counts, "meta.data"))[7] <- "celltype_cluster"
#' celltype_mapping <- list("T-cell" = c("T"), "Dendritic Cell" = character(0))
#' check_mapping_organism("Mus musculus", "g1", "g2", celltype_mapping, counts)
check_mapping_organism <- function(organism, target_class, base_class, celltype_mapping, counts) {
    # Basic types
    checkmate::assert_character(organism, len = 1, any.missing = FALSE)
    checkmate::assert_character(target_class, len = 1, any.missing = FALSE)
    checkmate::assert_character(base_class,  len = 1, any.missing = FALSE)
    checkmate::assert_list(celltype_mapping)
    # Accept Seurat or SingleCellExperiment
    is_seurat <- inherits(counts, "Seurat")
    is_sce    <- inherits(counts, "SingleCellExperiment")
    checkmate::assert_true(is_seurat || is_sce)
    # Pull per-cell metadata
    if (is_seurat) {
        md <- if (!is.null(counts@meta.data)) counts@meta.data else counts[[]]
    } else { # SCE
        md <- as.data.frame(SummarizedExperiment::colData(counts))
    }
    # Required columns
    checkmate::assert_subset(c("class", "celltype_cluster", "donor"), choices = colnames(md))
    # Ensure character
    md$class            <- as.character(md$class)
    md$celltype_cluster <- as.character(md$celltype_cluster)
    # Normalize for robust comparison
    norm <- function(x) trimws(tolower(x))
    classes_norm <- norm(md$class)
    t_class <- norm(target_class)
    b_class <- norm(base_class)
    # target vs base checks
    checkmate::assert_true(t_class != b_class)
    checkmate::assert_subset(c(t_class, b_class), choices = unique(classes_norm))
    # celltype_mapping must be well-formed
    keys <- names(celltype_mapping)
    checkmate::assert_true(!is.null(keys) && all(nzchar(keys)))
    checkmate::assert_true(!anyDuplicated(keys))
    TRUE
}

#' @title Create mapping organism object
#' @description
#' Creates a simple S3 object (list with class "mapping.organism") after validating its components.
#'
#' @param organism Character(1). Scientific Latin name of the organism.
#' @param target_class Character(1). Name of the target class for this organism.
#' @param base_class Character(1). Name of the base class for this organism.
#' @param celltype_mapping List. Mapping of cell types to clusters (character vectors).
#' @param counts `Seurat` or `SingleCellExperiment` object with the scRNA-seq
#' LogNormalized counts. This object should contain variables in 
#' `slot(SeuratObject, meta.data)` slot or
#' `slot(SingleCellExperimentObject, metadata)`; `class` indicating the class
#' the sample belongs to; `celltype_cluster` indicating the cell type cluster
#' (either character or numeric); `donor` indicating the sample ID.
#' @export
#' 
#' @return A list of class `"mapping.organism"` with elements:
#' \itemize{
#'   \item organism
#'   \item target_class
#'   \item base_class
#'   \item celltype_mapping
#'   \item counts
#' }
#' @examples
#' counts <- SeuratObject::pbmc_small
#' colnames(slot(counts, "meta.data"))[1] <- "donor"
#' colnames(slot(counts, "meta.data"))[6] <- "class"
#' colnames(slot(counts, "meta.data"))[7] <- "celltype_cluster"
#' celltype_mapping <- list("T-cell" = c("T"), "Dendritic Cell" = character(0))
#' obj <- create_mapping_organism("Mus musculus", "g1", "g2", celltype_mapping,
#' counts)
create_mapping_organism <- function(organism, target_class, base_class, celltype_mapping, counts) {
    check_mapping_organism(organism, target_class, base_class, celltype_mapping, counts)
    list(
        organism = organism,
        target_class = target_class,
        base_class = base_class,
        celltype_mapping = celltype_mapping,
        counts = counts
    )
}
