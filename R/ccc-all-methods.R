#' @title Pull Gene Set from MsigDB
#'
#' @description
#' \code{pullGeneSet} pulls the associated gene set from MsigDB. Default Gene
#' Set Collection (gse)cfrom MsigDB corresponds to homo sapiens (hs) genes with
#' ENTREZ ID and Human Gene Symbol (SYM) ID,cif any other gse is desired one
#' should load it with msigdb::getMsigdb(...).#'
#' @param object An object of class \link{pathway-class} or
#' \link{superpathway.gene.sets-class}. The object containing pathway
#' information to fetch gene sets.
#' @param gse  A Gene Set Collection (gse) from MsigDB. If `NULL` then gse
#' assigns homo sapiens organism with ENTREZ and HGNC IDs. It is recommended
#' to provide a gse object if multiple pathways are to be assessed, this will
#' reduce execution time.
#' @param ... Additional parameters passed to `msigdb::subsetCollection´ or
#' other related functions.
#'
#' @rdname pullGeneSet-method
#'
#' @return A vector containing the gene set corresponding to the specified
#' pathway
#' @import msigdb methods
#' @importFrom GSEABase geneIds
#'
#' @examples
#' # Example for pathway class
#' my_pathway <- new("pathway",
#'              standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#'              dbsource = "KEGG",
#'              collection = "c2",
#'              subcollection = "CP")
#' pullGeneSet(my_pathway)
#'
#' # Example for superpathway.gene.sets class
#' my_superpathway <- new("superpathway.gene.sets", pathway_info = my_pathway,
#'                        celltypes = c("T-cell", "Dendritic Cell")
#'                       )
#' pullGeneSet(my_superpathway)
methods::setGeneric("pullGeneSet",
                    function(object, gse = NULL, ...)
                        standardGeneric("pullGeneSet"))

#' @rdname pullGeneSet-method
#' @exportMethod pullGeneSet
#' @import msigdb methods
#' @importFrom GSEABase geneIds
#'
#' @examples
#' # Example for pathway class
#' my_pathway <- new("pathway",
#'              standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#'              dbsource = "KEGG",
#'              collection = "c2",
#'              subcollection = "CP")
#' pullGeneSet(my_pathway)
methods::setMethod("pullGeneSet",
                    "pathway",
                    function(object, gse = NULL, ...){
                        if(is.null(gse)){
                            gse <- msigdb::getMsigdb(
                                org = c("hs"), id = c("SYM", "EZID"), 
                                version = msigdb::getMsigdbVersions()[1])
                            }
                        # Database to grab gene set from
                        database <- paste0(
                            object@subcollection, ":", object@dbsource)
                        # If database is KEGG one should append KEGG
                        if(object@dbsource == "KEGG"){
                            gse_collection_full <- msigdb::subsetCollection(
                                gse, collection = object@collection,
                                subcollection = object@subcollection,  ...)
                            gse_collection_append <- msigdb::appendKEGG(
                                gse_collection_full,
                                version = msigdb::getMsigdbVersions())
                            gse_collection_kegg <- msigdb::subsetCollection(
                                gse_collection_append,
                                collection = object@collection,
                                subcollection = database)
                            Pathway_info <-
                                gse_collection_kegg[[object@standard_name]]
                            Pathway_gene_set <- GSEABase::geneIds(Pathway_info)
                            }else{
                                gse_collection_full <- msigdb::subsetCollection(
                                    gse, collection = object@collection,
                                    subcollection = database, ...)
                                Pathway_info <- gse_collection_full[[
                                    object@standard_name]]
                                Pathway_gene_set <- GSEABase::geneIds(
                                    Pathway_info)
                                }
                        return(Pathway_gene_set)
                        }
)

#' @rdname pullGeneSet-method
#' @exportMethod pullGeneSet
#'
#' @import msigdb methods
#' @importFrom GSEABase geneIds
#'
#' @examples
#' # Example for superpathway.gene.sets class
#' my_pathway <- new("pathway",
#' standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#' dbsource = "KEGG",
#' collection = "c2",
#' subcollection = "CP")
#' my_superpathway <- new("superpathway.gene.sets", pathway_info = my_pathway,
#'                        celltypes = c("T-cell", "Dendritic Cell")
#'                        )
#' pullGeneSet(my_superpathway)
methods::setMethod("pullGeneSet",
                    "superpathway.gene.sets",
                    function(object, gse = NULL, ...){
                        if(is.null(gse)){
                            gse <- msigdb::getMsigdb(
                                org = c("hs"), id = c("SYM", "EZID"), 
                                version = msigdb::getMsigdbVersions()[1])
                            }
                        # Database to grab gene set from
                        database <- paste0(
                            object@pathway_info@subcollection, ":",
                            object@pathway_info@dbsource)
                        # If database is KEGG one should append KEGG
                        if(object@pathway_info@dbsource == "KEGG"){
                            gse_collection_full <- msigdb::subsetCollection(
                                gse,
                                collection = object@pathway_info@collection,
                                subcollection=object@pathway_info@subcollection,
                                ...)
                            gse_collection_append <- msigdb::appendKEGG(
                                gse_collection_full,
                                version = msigdb::getMsigdbVersions())
                            gse_collection_kegg <- msigdb::subsetCollection(
                                gse_collection_append,
                                collection = object@pathway_info@collection,
                                subcollection = database)
                            Pathway_info <- gse_collection_kegg[[
                                object@pathway_info@standard_name]]
                            Pathway_gene_set <- GSEABase::geneIds(Pathway_info)
                            }else{
                                gse_collection_full <- msigdb::subsetCollection(
                                    gse,
                                    collection = object@pathway_info@collection,
                                    subcollection = database, ...)
                                Pathway_info <- gse_collection_full[[
                                    object@pathway_info@standard_name]]
                                Pathway_gene_set<-
                                    GSEABase::geneIds(Pathway_info)
                                }
                        return(Pathway_gene_set)
                        }
)

#' @title Setter for gene_sets_celltype of superpathway.gene.sets Class
#'
#' @description
#' A setter for gene_sets_celltype slot of superpathway.gene.sets Class that
#' checks for its validity when updating. Number of gene sets should be equal to
#'  the number of cell types, the updated slot cannot have NULL value
#'
#'
#' @param x A \link{superpathway.gene.sets-class} object to update
#' @param value A list with the gene sets per cell type
#'
#' @exportMethod setGeneSetsCelltype<-
#'
#' @rdname setGeneSetsCelltype-method
#' @return A \link{superpathway.gene.sets-class} object with the updated
#' gene_sets_celltype slot
methods::setGeneric("setGeneSetsCelltype<-",
                    function(x, value) standardGeneric("setGeneSetsCelltype<-"))

#' @rdname setGeneSetsCelltype-method
#' @exportMethod setGeneSetsCelltype<-
#'
#' @import msigdb methods
#' @importFrom GSEABase geneIds
#'
#' @examples
#' my_pathway <- new("pathway",
#'                standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#'                 dbsource = "KEGG",
#'                 collection = "c2",
#'                 subcollection = "CP")
#' # Superpatheway.gene.sets with NULL gene_sets_celltype
#' my_superpathway <- new("superpathway.gene.sets", pathway_info = my_pathway,
#'                        celltypes = c("T-cell", "Dendritic Cell")
#'                        )
#' setGeneSetsCelltype(my_superpathway) <- list(c("IL2"), c("IL4"))
methods::setMethod("setGeneSetsCelltype<-",
                    "superpathway.gene.sets",
                    function(x, value){
                    # Assess validity conditions of value
                    checkmate::assert_list(value, null.ok = FALSE)
                    if(!is.null(value)){
                        len_genesets <- length(value)
                        len_celltypes <- length(x@celltypes)
                        checkmate::assert_true(len_genesets == len_celltypes)
                        }
                        x@gene_sets_celltype <- value
                        return(x)
                        }
)

#' @title Repeat gene sets per cell type
#'
#' @description
#' \code{setRepeatGeneSets} pulls the gene set from MsigDB via
#' \link{pullGeneSet} and repeats it as many times as cell types. This method is
#'  useful when all cell types have the same gene set.
#' @param object An object of class \link{superpathway.gene.sets-class} to
#' assign the repeated gene sets per cell type.
#' @param ... Other parameters to pass onto \link{pullGeneSet}
#' 
#' @rdname setRepeatGeneSets-method
#'
#' @return A \link{superpathway.gene.sets-class} object with updated
#' gene_sets_celltype slot with the repeated gene sets.
methods::setGeneric("setRepeatGeneSets",
                    function(object, ...) standardGeneric("setRepeatGeneSets"))

#' @rdname setRepeatGeneSets-method
#' @exportMethod setRepeatGeneSets
#' @examples
#' my_pathway <- new("pathway",
#'              standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#'              dbsource = "KEGG",
#'              collection = "c2",
#'              subcollection = "CP")
#' my_superpathway <- new("superpathway.gene.sets", pathway_info = my_pathway,
#'                        celltypes = c("T-cell", "Dendritic Cell")
#'                        )
#' my_superpathway <- setRepeatGeneSets(my_superpathway)
#' print(my_superpathway)
methods::setMethod("setRepeatGeneSets",
                    "superpathway.gene.sets",
                    function(object, ...){
                        # Pull Gene Set
                        pulled_gene_set <- pullGeneSet(object, ...)
                        repeated_gene_set <- rep(list(pulled_gene_set),
                                            length(object@celltypes))
                        # Update superpathway.gene.sets object
                        setGeneSetsCelltype(object) <- repeated_gene_set
                        return(object)
                    }
)

#' @title Build predictor block and response matrix
#'
#' @description
#' \code{matrixToBlock} generates the predicted block with the
#' pseudobulk_lognorm matrix, the response matrix with the sample_class slot and
#' the dimensions of each block, with the information provided in a
#' superpathway.input object.
#'
#' @param object An object of class \link{superpathway.input-class}.
#'
#' @rdname matrixToBlock-method
#'
#' @return A list containing the predictor block, response matrix and block
#' dimension for fit asmbPLS-DA
#' @exportMethod matrixToBlock
#' @import methods
methods::setGeneric("matrixToBlock",
                    function(object) standardGeneric("matrixToBlock"))

#' @rdname matrixToBlock-method
#' @importFrom methods signature
#' @exportMethod matrixToBlock
methods::setMethod("matrixToBlock",
                    signature = methods::signature(
                        object = "superpathway.input"),
                    definition = matrixToBlock.superpathway.input)

#' @title Fit optimal asmbPLS-DA
#'
#' @param object A \link{superpathway.input-class} object
#' @param ... Other parameters to pass to \link{fitOptimal.superpathway.input}
#'
#' @rdname fitOptimal-method
#'
#' @exportMethod fitOptimal
methods::setGeneric("fitOptimal",
                    function(object, ...) standardGeneric("fitOptimal"))

#' @rdname fitOptimal-method
#' @importFrom methods signature
#' @exportMethod fitOptimal
methods::setMethod("fitOptimal",
                    signature = methods::signature(
                        object = "superpathway.input"),
                    definition = fitOptimal.superpathway.input)

