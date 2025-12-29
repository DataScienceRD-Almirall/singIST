#' @title Pull Gene Set from MsigDB
#' @description
#' Retrieves the gene set associated with a pathway or superpathway from MsigDB.
#'
#' @param object List. A pathway (from create_pathway()) or superpathway (from create_superpathway()).
#' @param gse Gene Set Collection from MsigDB. If NULL, loads default (human, ENTREZ + SYM IDs).
#' @param ... Additional arguments passed to msigdb::subsetCollection().
#'
#' @return A character vector with gene IDs for the specified pathway.
#' @import msigdb
#' @importFrom GSEABase geneIds
#' @export
#' @examples
#' \donttest{
#' my_pathway <- create_pathway("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#' "KEGG", "c2", "CP")
#' pullGeneSet(my_pathway)}
pullGeneSet <- function(object, gse = NULL, ...) {
    if (is.null(gse)) {
        gse <- msigdb::getMsigdb(org = "hs", id = c("SYM", "EZID"),
                                 version = msigdb::getMsigdbVersions()[1])
    }
    # Detect if object is pathway or superpathway
    if (!is.null(object$pathway_info)) {
        pw <- object$pathway_info
    } else {
        pw <- object
    }
    database <- paste0(pw$subcollection, ":", pw$dbsource)
    if (pw$dbsource == "KEGG") {
        gse_full <- msigdb::subsetCollection(gse, collection = pw$collection,
                                             subcollection = pw$subcollection, ...)
        gse_kegg <- msigdb::subsetCollection(msigdb::appendKEGG(gse_full,
                                                                version = msigdb::getMsigdbVersions()),
                                             collection = pw$collection,
                                             subcollection = database)
        pathway_info <- gse_kegg[[pw$standard_name]]
    } else {
        gse_full <- msigdb::subsetCollection(gse, collection = pw$collection,
                                             subcollection = database, ...)
        pathway_info <- gse_full[[pw$standard_name]]
    }
    return(GSEABase::geneIds(pathway_info))
}

#' @title Set gene sets per cell type in a superpathway
#' @description
#' Updates the gene_sets_celltype element of a superpathway object,
#' ensuring validity:
#' - The number of gene sets must match the number of cell types.
#'
#' @param object List. A superpathway object created by create_superpathway(), or
#' a pathway object created by create_pathway().
#' @param value List. A list of genes to incorporate in each cell type slot. By
#' default `NULL`, only use if genes are to be introduced manually.
#' @param ... other parameters to pass onto pullGeneSet().
#' 
#' @return The updated superpathway object (list).
#' @export
#' @examples
#' my_pathway <- create_pathway(
#'   standard_name = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#'   dbsource = "KEGG",
#'   collection = "c2",
#'   subcollection = "CP"
#' )
#' 
#' my_superpathway <- create_superpathway(my_pathway, c("T-cell",
#' "Dendritic Cell"), list())
#' \donttest{
#' my_superpathway <- setGeneSetsCelltype(my_superpathway, list(c("IL2"),
#' c("IL4")))}
setGeneSetsCelltype <- function(object, value = NULL, ...) {
    # Check if list is a pathway or superpathway type
    if("gene_sets_celltype" %in% names(object)){
        len_celltypes <- length(object$celltypes)
        if(is.null(value)){
            set <- pullGeneSet(object$pathway_info, ...)
            object$gene_sets_celltype <- rep(list(set), len_celltypes)
        }else{
            set <- value
        }
    }else{
        stop("This is not a superpathway object")
    }
    return(object)
}
