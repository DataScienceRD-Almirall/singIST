#' @title Pathway information
#'
#' @description
#' Pathway from MsigDB to grab gene set from.
#'
#' @slot standard_name A character with the pathway standard name provided by
#' MsigDB.
#' @slot dbsource A character with the pathway database to grab information
#' from. Currently available options are: KEGG, PID, REACTOME, BIOCARTA,
#' WP (WikiPathways).
#' @slot collection A character with MsigDB collection to grab information from
#' (currently only holds c2).
#' @slot subcollection A character with MsigDB subcollection to grab information
#' from (currently only holds CP).
#'
#' @name pathway-class
#' @rdname pathway-class
#'
#' @exportClass pathway
#'
#' @import methods checkmate
#'
#' @examples
#' my_pathway <- new("pathway",
#' standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#' dbsource = "KEGG",
#' collection = "c2",
#' subcollection = "CP")
#' print(my_pathway)
methods::setClass("pathway",
                    slots = list(standard_name = "character",
                                dbsource = "character",
                                collection = "character",
                                subcollection = "character"
                            ),
                    validity = function(object){
                        checkmate::assert_character(object@standard_name)
                        # Currently only KEGG, PID, REACTOME, BIOCARTA and WP
                        # pathway databases
                        checkmate::assert_choice(
                                    object@dbsource,
                                    choices = c("KEGG", "PID", "REACTOME",
                                                "BIOCARTA", "WP"))
                    # Currently only Curated Gene sets (C2) from Canonical
                    # Pathways (CP) are possible
                    checkmate::assert_choice(
                        object@collection,
                        choices = c("c2"))
                    checkmate::assert_choice(
                        object@subcollection,
                        choices = c("CP"))
                    # Check that standard name is consistent with the provided
                    # dbsource
                    checkmate::assert_character(
                        object@standard_name,
                        pattern = paste0("^", object@dbsource, "_"))
                    return(TRUE)
                    }
)

#' @title Superpathway gene sets per cell type
#'
#' @description
#' Class containing the superpathway information: cell types and gene sets per
#' each cell type.
#'
#' @slot pathway_info A \link{pathway-class} object.
#' @slot celltypes A vector where each element is a character representing the
#' cell type.
#' @slot gene_sets_celltype A list of character vectors, each list element
#' corresponds to a gene set for each cell type. If the gene sets are identical
#' for all cell types look at method MIRACOÑO.
#'
#' @name superpathway.gene.sets-class
#' @rdname superpathway.gene.sets-class
#'
#' @import methods checkmate
#' @exportClass superpathway.gene.sets
#'
#' @examples
#' my_pathway <- new("pathway",
#' standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#' dbsource = "KEGG",
#' collection = "c2",
#' subcollection = "CP")
#'
#' my_superpathway <- new("superpathway.gene.sets",
#'                  pathway_info = my_pathway,
#'                  celltypes = c("T-cell", "Dendritic Cell"),
#'                  gene_sets_celltype = list(c("IL4", "IL5"), c("IL13")))
#' print(my_superpathway)
methods::setClass("superpathway.gene.sets",
                    slots = list(pathway_info = "pathway",
                                celltypes = "character",
                                gene_sets_celltype = "list"),
                    validity = function(object) {
                        checkmate::assert_class(
                            object@pathway_info,
                            "pathway")
                        checkmate::assert_character(
                            object@celltypes,
                            min.len = 2)
                        checkmate::assert_list(
                            object@gene_sets_celltype,
                            null.ok = TRUE)
                        if(!is.null(object@gene_sets_celltype)){
                            len_genesets <- length(object@gene_sets_celltype)
                            len_celltypes <- length(object@gene_sets_celltype)
                            checkmate::assert_true(len_genesets==len_celltypes)
                        }
                        return(TRUE)
                    }
)

#' @title asmbPLS-DA hyperparameters
#'
#' @description
#' A class with the hyperparameters for fit and cross-validation procedure of
#' asmbPLS-DA.
#'
#' @slot quantile_comb_table A matrix containing the quantile (lambda) sparsity
#' values for the asmbPLS-DA cross validation step. Rows define combination of
#' quantiles and columns define cell types.
#' @slot outcome_type A character indicating whether output is "binary" or
#' "multiclass" for asmbPLS-DA
#' @slot number_PLS An integer indicating the maximum number of PLS components
#' for asmbPLS-DA
#' @slot folds_CV An integer indicating the number of folds for Cross
#' Validation. If NULL, the default is 5. If the number is 1 a Leave-One-Out
#' Cross Validation (LOOCV) is performed. LOOCV is automatically performed if
#' the number of samples per class is less than 5, indistinguisably of folds_CV
#' number.
#' @slot repetition_CV An integer indicating the number of repetitions of CV.
#' If NULL, the default is 10. If folds_CV is 1 then repetition_CV should also
#' be 1.
#'
#' @name hyperparameters-class
#' @rdname hyperparameters-class
#'
#' @import methods checkmate
#' @importFrom RcppAlgos permuteGeneral
#' @exportClass hyperparameters
#' @examples
#' quantile_comb_table <- base::as.matrix(
#' RcppAlgos::permuteGeneral(seq(0.05, 0.95, by = 0.50)),
#' ncol = length(c("T-cell", "Dendritic Cell"))
#' )
#' outcome_type <- "binary"
#' number_PLS <- as.integer(3)
#' folds_CV <- as.integer(1)
#' repetition_CV <- as.integer(1)
#' my_hyperparameters <- new("hyperparameters",
#'                            quantile_comb_table = quantile_comb_table,
#'                            outcome_type = outcome_type,
#'                            number_PLS = number_PLS,
#'                            folds_CV = folds_CV
#'                            )
#' print(my_hyperparameters)
methods::setClass("hyperparameters",
                    slots = list(
                    quantile_comb_table = "matrix",
                    outcome_type = "character",
                    number_PLS = "integer",
                    folds_CV = "integer",
                    repetition_CV = "integer"
                    ),
                    validity = function(object){
                        checkmate::assert_matrix(
                            object@quantile_comb_table, min.rows = 1)
                        checkmate::assert_choice(
                            object@outcome_type,
                            choices = c("binary", "multiclass"))
                        checkmate::assert_integer(
                            object@number_PLS, lower = 0)
                        checkmate::assert_integer(
                            object@folds_CV, lower = 0, null.ok = TRUE)
                        checkmate::assert_integer(
                            object@repetition_CV, lower = 0, null.ok = TRUE)
                        return(TRUE)
                    }
)

#' @title Superpathway input for asmbPLS-DA
#'
#' @slot superpathway_info A \link{superpathway.gene.sets-class} object
#' @slot hyperparameters_info A \link{hyperparameters-class} object
#' @slot pseudobulk_lognorm A matrix from Seurat::AggregateExpression() where
#' columns are genes and rows are combinations of sample id and cell type in
#' the form "Celltype_Sampleid". Rownames should be of the form
#' "Celltype_Sampleid" and columns of the form of "HGNC".
#' @slot sample_id A vector of characters with the sample id
#' @slot sample_class A vector of characters with the class of each sample
#' @slot base_class A character indicating the base class
#' @slot target_class A character indicating the target class
#'
#' @name superpathway.input-class
#' @rdname superpathway.input-class
#'
#' @import methods checkmate
#' @importFrom RcppAlgos permuteGeneral
#' @exportClass superpathway.input
#'
#' @examples
#' # Load pathway and superpathway.gene.sets objects
#' my_pathway <- new("pathway",
#' standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#' dbsource = "KEGG",
#' collection = "c2",
#' subcollection = "CP")
#'
#' celltypes <- c("T-cell", "Dendritic Cell")
#' my_superpathway <- new("superpathway.gene.sets",
#'                  pathway_info = my_pathway,
#'                  celltypes = celltypes)
#'
#' # Load hyperparameters object
#' quantile_comb_table <- base::as.matrix(
#' RcppAlgos::permuteGeneral(seq(0.05, 0.95, by = 0.50)),
#'ncol = length(celltypes)
#')
#' outcome_type <- "binary"
#' number_PLS <- as.integer(3)
#' folds_CV <- as.integer(1)
#' repetition_CV <- as.integer(1)
#' my_hyperparameters <- new("hyperparameters",
#'                            quantile_comb_table = quantile_comb_table,
#'                            outcome_type = outcome_type,
#'                            number_PLS = number_PLS,
#'                            folds_CV = folds_CV
#'                            )
#' # Slots of superpathway input
#' sample_id <- c("AD1", "AD2", "HC1", "HC2")
#' sample_class <- c("AD", "AD", "HC", "HC")
#' base_class <- "HC"
#' target_class <- "AD"
#' pseudobulk_lognorm <- matrix(rnorm(length(celltypes)*length(sample_id)),
#' nrow = length(celltypes)*length(sample_id), ncol = length(celltypes))
#' rownames(pseudobulk_lognorm) <- as.vector(t(outer(celltypes, sample_id,
#' function(x, y) paste(x, y, sep = "_"))))
#'
#' my_superpathway_input <- new("superpathway.input",
#' superpathway_info = my_superpathway,
#' hyperparameters_info = my_hyperparameters,
#' pseudobulk_lognorm = pseudobulk_lognorm,
#' sample_id = sample_id,
#' sample_class = sample_class,
#' base_class = base_class,
#' target_class = target_class)
#' print(my_superpathway_input)
methods::setClass("superpathway.input",
                    slots = list(
                        superpathway_info = "superpathway.gene.sets",
                        hyperparameters_info = "hyperparameters",
                        pseudobulk_lognorm = "matrix",
                        sample_id = "character",
                        sample_class = "character",
                        base_class = "character",
                        target_class = "character"
                    ),
                    validity = function(object){
                        checkmate::assert_class(
                            object@superpathway_info, "superpathway.gene.sets")
                    checkmate::assert_class(
                        object@hyperparameters_info, "hyperparameters")
                    checkmate::assert_character(object@sample_id)
                    checkmate::assert_character(object@sample_class)
                    checkmate::assert_choice(
                        object@base_class,
                        choices = unique(object@sample_class))
                    possible_classes <- setdiff(
                        unique(object@sample_class), object@base_class)
                    checkmate::assert_choice(
                        object@target_class, choices = possible_classes)
                    # Check for consistency of the pseudobulk matrix with
                    # sample_id and celltypes information
                    checkmate::assert_matrix(
                        object@pseudobulk_lognorm, any.missing = FALSE)
                    checkmate::assert_true(all(as.vector(
                        t(outer(object@superpathway_info@celltypes,
                                object@sample_id,
                                function(x, y)
                                    paste(x, y, sep = "_")))) ==
                                    rownames(object@pseudobulk_lognorm)))
                    # Check that gene_sets is not NULL
                    checkmate::assert_list(
                        object@superpathway_info@gene_sets_celltype)
                    # Check consistency between
                    # hyperparameters_info@outcome.type with the
                    # provided sample_class
                    if(object@hyperparameters_info@outcome_type == "binary"){
                        checkmate::assert_true(
                            length(unique(object@sample_class)) == 2)
                    }else{
                        checkmate::assert_true(
                            length(unique(object@sample_class)) >= 3)
                    }
                    # Check that sample id does not contain "_" characters
                    sample_character <- base::grepl("_", object@sample_id)
                    checkmate::assert_true(!all(sample_character))
                    return(TRUE)
                    }
)

#' @title Fit asmbPLS-DA with superpathway input
#'
#' @slot superpathway_input A \link{superpathway.input-class} object with the
#' superpathway input used to fit asmbPLS-DA
#' @slot hyperparameters_fit A \link{hyperparameters-class} object with the
#' hyperparameters used to fit
#' @slot model_fit A list with the fitted model
#' @slot model_validation A list with the validation metrics of the fitted model
#'
#' @name superpathway.fit.model-class
#' @rdname superpathway.fit.model-class
#'
#' @import methods checkmate
#' @exportClass superpathway.fit.model
methods::setClass("superpathway.fit.model",
                    slots = list(
                    superpathway_input = "superpathway.input",
                    hyperparameters_fit = "hyperparameters",
                    model_fit = "list",
                    model_validation = "list"
                    ),
                validity = function(object){
                    checkmate::assert_class(
                        object@superpathway_input, "superpathway.input")
                    checkmate::assert_class(
                        object@hyperparameters_fit, "hyperparameters")
                    checkmate::assert_list(
                        object@model_fit)
                    checkmate::assert_list(
                        object@model_validation)
                    return(TRUE)
                    }
)

#' @title Mapping organism input
#'
#' @slot organism A character with the scientific Latin name of the organism
#' @slot target_class A character indicating the name of the target class for
#' this organism
#' @slot base_class A character indicating the name of the base class for this
#' organism
#' @slot celltype_mapping A list of vectors with the cell type correspondence
#' between the mapping organism and the reference organism for which asmbPLSDA
#' has been trained. Note that the name of each element of the list should be
#' exactly the cell type used in \link{superpathway.gene.sets-class}, each
#' vector contains the values of `slot(SeuratObject, meta.data$class)` that are
#' equivalent to its respective \link{superpathway.gene.sets-class} cell type.
#' If no mapping exists for a given cell type its vector should void. If you
#' are assessing multiple \link{superpathway.gene.sets-class} objects, you
#' should include the mapping of all cell types used in these objects.
#' @slot counts A Seurat object with the scRNA-seq counts. This object should
#' contain variables in `slot(SeuratObject, meta.data)` slot; `class`
#' indicating the class the sample belongs to; `celltype_cluster` indicating
#' the cell type cluster (either character or numeric)
#'
#' @name mapping.organism-class
#' @rdname mapping.organism-class
#'
#' @exportClass mapping.organism
#'
#' @import checkmate SeuratObject sp
#' @examples
#' organism <- "Mus musculus"
#' target_class <- "g1"
#' base_class <- "g2"
#' counts <- SeuratObject::pbmc_small # Toy dataset
#' # Rename "group" variable to "class"
#' colnames(slot(counts, "meta.data"))[6] <- "class"
#' # Example existing mapping for T-cell but no mapping for Dendritic Cell
#' celltype_mapping <- list("T-cell" = c("dγdT", "T"),
#' "Dendritic Cell" = c())
#' # Rename "RNA_snn.res.1" variable to "celltype_cluster"
#' colnames(slot(counts, "meta.data"))[7] <- "celltype_cluster"
#' # Create object
#' my_mapping_organism <- new("mapping.organism",
#'                            organism = organism,
#'                            target_class = target_class,
#'                            base_class = base_class,
#'                            celltype_mapping = celltype_mapping,
#'                            counts = counts)
#' print(my_mapping_organism)
methods::setClass("mapping.organism",
                    slots = list(
                        organism = "character",
                        target_class = "character",
                        base_class = "character",
                        celltype_mapping = "list",
                        counts = "Seurat"
                    ),
                    validity = function(object){
                        checkmate::assert_class(object@counts, "Seurat")
                        checkmate::assert_character(object@organism)
                        checkmate::assert_list(object@celltype_mapping)
                        # Check that class and celltype_cluster columns exist
                        checkmate::assert_true(
                            all(c("class", "celltype_cluster") %in%
                                    colnames(object@counts@meta.data)))
                        # Check that target and base classes exists in the
                        # variable class
                        checkmate::assert_true(
                            all(c(object@target_class, object@base_class) %in%
                                    object@counts@meta.data$class))
                        # Check that target class is different from base class
                        checkmate::assert_true(
                            object@target_class != object@base_class)
                        return(TRUE)
                    }
)
