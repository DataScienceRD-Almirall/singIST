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
#' file <- system.file("extdata", "example_mapping_organism.rda",
#' package = "singIST")
#' load(file)
#' data <- example_mapping_organism
#' new_object <- celltype_mapping(data)
#' slot(new_object, "counts")$celltype_cluster
celltype_mapping <- function(object){
    checkmate::assert_class(object, "mapping.organism")
    output <- object
    # Avoid spaces as FindMarkers/findMarkers does not identify them
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

#' @title Compute differentially expressed genes for pseudobulk with limma-voom
#' @description
#' Computes differentially expressed genes with limma-voom after pseudobulking
#' the scRNAseq matrix. The reported logFC values are difference of means of
#' log-normalized expression values with `Seurat::AggregateExpression` or
#' `SingleCellExperiment::aggregateAcrossCells`. This logFC is consistent with
#' the human log2FC computation by asmbPLS-DA. 
#' @param object A \link{mapping.organism-class} object. 
#' @param logfc.treshold Sets the minimum log-fold change (logFC) cutoff for
#' identifying differentially expressed genes (DEGs). By default
#' `logfc.treshold = 0`.
#' @param ... Other parameters to pass onto `limma::voom` or `limma::lmFit`.
#' @param assay Specific assay being used for analysis. By default
#' `assay = RNA`.
#' @import checkmate Seurat scran SingleCellExperiment
#' @returns
#' A list where each element is a data.frame for a cell type containing: `p_val`
#' p-value of limma-voom test, `avg_log2FC` descriptive point estimate of
#' logFC, and `p_val_adj` FDR for all the transcriptome as output by limma-voom. 
#' @import Seurat edgeR limma dplyr
#' @export
#' @examples
#' # Set the identities
#' file <- system.file("extdata", "example_mapping_organism.rda",
#' package = "singIST")
#' load(file)
#' data_organism <- example_mapping_organism
#' data <- celltype_mapping(data_organism)
#' diff_expressed(data)
diff_expressed <- function(object, assay = "RNA", logfc.treshold = 0, ...){
    checkmate::assert_class(object, "mapping.organism")
    # 1) PSEUDOBULK -----------------------------------------------------------
    sce <- object@counts
    if (inherits(sce, "Seurat")) {
        pb <- Seurat::AggregateExpression(
            object = sce,
            assays = assay,
            group.by = c("celltype_cluster", "class", "donor"),
            return.seurat = TRUE
        )
        mat <- round(Seurat::GetAssayData(pb, assay=assay, slot="counts"))
        meta <- pb@meta.data[, c("celltype_cluster", "class", "donor"),
                             drop=FALSE]
    } else if (inherits(sce, "SingleCellExperiment")) {
        require(scuttle)
        grp <- paste0(
            SummarizedExperiment::colData(sce)$celltype_cluster, "_",
            SummarizedExperiment::colData(sce)$class, "_",
            SummarizedExperiment::colData(sce)$donor
        )
        pb_sce <- scuttle::aggregateAcrossCells(sce, ids=DataFrame(grp=grp))
        mat <- assay(pb_sce, assay)
        tmp <- strsplit(as.character(colData(pb_sce)$grp), "_", fixed=TRUE)
        meta <- data.frame(
            celltype_cluster = vapply(tmp, `[`, 1, FUN.VALUE=""),
            class = vapply(tmp, `[`, 2, FUN.VALUE=""),
            donor = vapply(tmp, `[`, 3, FUN.VALUE=""),
            row.names = colnames(mat),
            stringsAsFactors = FALSE
        )
    }
    # 2) BUILD LIMMA–VOOM MODEL ACROSS ALL PSEUDOBULK SAMPLES -----------------
    meta$celltype_cluster <- gsub("-", "_",meta$celltype_cluster)
    meta$class <- gsub("-", "_", meta$class)
    meta$grp <- factor(paste(meta$celltype_cluster, meta$class, sep="_"))
    dge <- edgeR::DGEList(counts = mat)
    dge <- edgeR::calcNormFactors(dge)
    design <- stats::model.matrix(~0 + grp, data=meta)
    colnames(design) <- levels(meta$grp)
    v <- limma::voom(dge, design, plot=FALSE, ...)
    fit <- limma::lmFit(v, design, ...)
    # for each cell-type, define the "target-base" contrast
    celltypes <- sort(unique(meta$celltype_cluster))
    contrasts <- paste0(
        celltypes, "_", object@target_class, 
        " - ",
        celltypes, "_", object@base_class
    )
    cm <- makeContrasts(contrasts=contrasts, levels=design)
    fit2 <- contrasts.fit(fit, cm) %>% limma::eBayes()
    # 3) EXTRACT FOR PER CELL TYPE ------------------------------
    all_pvals <- fit2$p.value # genes × contrasts
    genes <- rownames(fit2)
    out <- setNames(vector("list", length(celltypes)), celltypes)
    for (i in seq_along(celltypes)) {
        b <- celltypes[i]
        cn <- colnames(cm)[i] # e.g. "Bcell_target - Bcell_base"
        pv <- all_pvals[genes, cn]
        pv_adj <- p.adjust(pv, method="BH") # 
        # descriptive log2FC on voom E (log2-CPM)
        idx_t <- meta$celltype_cluster==b & meta$class==object@target_class
        idx_b <- meta$celltype_cluster==b & meta$class==object@base_class
        if(inherits(pb, "Seurat")){
            mat_log <- Seurat::GetAssayData(pb[genes, ],
                                            assay=assay,slot="data")
        }
        descr_fc <- rowMeans(mat_log[genes, idx_t, drop = FALSE]) -
            rowMeans(mat_log[genes, idx_b, drop = FALSE])
        df <- data.frame(
            "p_val" = pv,
            "avg_log2FC" = descr_fc,
            "pct.1" = rep("Not available", length(descr_fc)),
            "pct.2" = rep("Not available", length(descr_fc)),
            "p_val_adj" = pv_adj, # This will be corrected only with gene set
            row.names = genes,
            stringsAsFactors = FALSE
        )
        out[[b]] <- df
    }
    return(out)
}

#' title Compute differentially expressed genes with FindMarkers/findMarkers
#' description
#' Computes differentially expressed genes with Seurat::FindMarkers, for
#' Seurat objects, or scran::findMarkers, for SingleCellExperiment
#' objects, for the conditions indicated. Note that Seurat::FindMarkers will
#' compute Wilcoxon Signed Rank Test by default, while scran::findMarkers will
#' perform t-test by default instead. The reported logFC values are differences
#' in the log of the average exponentiated data, with pseudocount. This logFC
#' reporting is consistent with Seurat::FindMarkers, and differs from methods
#' like scran::findMarkers, which compute logFC as the difference in mean
#' log expression values between groups. Hence, logFC will be differently
#' reported than scran::findMarkers to ensure reproducibility and consistency
#' independently of the class, Seurat or SingleCellExperiment, provided. The
#' output is rendered to be homogeneous indistinguishably of the method used.
#' 
#' param object A link{mapping.organism-class} object. If a Seurat object
#' was provided, then Idents(object) assigned to variables with the conditions
#' being tested is expected.
#' param condition_1 A vector with the elements of the first factor to perform
#' the hypothesis test. By default the mapped cell types
#' condition_1 = names(slot(object, "celltype_mapping"))
#' param condition_2 A vector with the elements of the second factor to perform
#' the hypothesis test with. By default the class of the organism
#' condition_2 = c(slot(object, "target_class"), slot(object, "base_class"))
#' param logfc.treshold Sets the minimum log-fold change (logFC) cutoff for
#' identifying differentially expressed genes (DEGs). By default
#' logfc.treshold = 0.25
#' param ... Other parameters to pass onto Seurat::FindMarkers() or
#' scran::findMarkers.
#' param assay Specific assay being used for analysis. By default
#' assay = RNA.
#' import checkmate Seurat scran SingleCellExperiment
#' returns
#' A list where each element is a data.frame for a cell type containing: p_val
#' p-value of test, avg_log2FC logFC between mean expression values,
#' pct.1 percentage of cells where the gene is detected in the base class,
#' pct.2 percentage of cells where the gene is detected in the target class,
#' p_val_adj FDR.
#' export
#' examples
#' # Set the identities
#' file <- system.file("extdata", "example_mapping_organism.rda",
#' package = "singIST")
#' load(file)
#' data_organism <- example_mapping_organism
#' data <- celltype_mapping(data_organism)
#' slot(data, "counts")$test <- paste0(slot(data, "counts")$celltype_cluster,
#' "_", slot(data, "counts")$class)
#' SeuratObject::Idents(slot(data, "counts")) <- "test"
#' diff_expressed(data)
# diff_expressed <- function(object, condition_1 = c(), condition_2 = c(),
#                             logfc.treshold = 0.25, assay = "RNA", ...){
#     checkmate::assert_class(object, "mapping.organism")
#     if(is.null(condition_1)){condition_1 <-
#         names(object@celltype_mapping)[lengths(object@celltype_mapping)>0]}
#     if(is.null(condition_2)){condition_2 <-
#         c(object@target_class, object@base_class)}
#     counts <- object@counts
#     if(is(object@counts,"Seurat")){
#         counts_aggr <- AggregateExpression(counts, assays = "RNA", slot = "data",
#                                            return.seurat = TRUE,
#                                            group.by = c("celltype_cluster",
#                                                         "class",
#                                                         "donor"))
#         counts_aggr$celltype_cluster <- gsub("-", "_", counts_aggr$celltype_cluster)
#         counts_aggr$test <- paste0(counts_aggr$celltype_cluster, "_", counts_aggr$class)
#         SeuratObject::Idents(counts_aggr) <- "test"
#         #counts_aggr[["RNA"]]$counts <- round(counts_aggr[["RNA"]]$counts)
#         apply_function <- function(row, data = counts_aggr, ...) {
#             logFC <- Seurat::FindMarkers(
#                 object = data, ident.1 = row[1], ident.2 = row[2], slot = "data",
#                 logfc.threshold = logfc.treshold,
#                 test.use = "t", ...)
#             return(logFC)
#         }
#         # Combinations to test
#         combinations <- base::outer(condition_1, condition_2 , paste, sep = "_")
#         output <- base::apply(combinations, 1, apply_function)
#         names(output) <- condition_1
#     }else{ # If not Seurat then its SingleCellExperiment
#         output <-lapply(condition_1, function(x, sce=counts, lfc=logfc.treshold,
#                                                 classes = condition_2, ...){
#             filt <- sce[, sce$celltype_cluster == x]
#             DEG <- scran::findMarkers(filt, groups = filt$class, lfc = lfc,
#                                         add.summary = TRUE, min.prop = 0.01,
#                                         assay.type = "logcounts", ...)
#             norm_expr <- expm1(SingleCellExperiment::logcounts(filt))
#             group_base <- colnames(filt)[filt$class == condition_2[2]]
#             group_target <- colnames(filt)[filt$class == condition_2[1]]
#             mean_base <- log(
#                 (rowSums(norm_expr[, group_base]) + 1)/length(group_base),
#                 base = 2)
#             mean_target <- log(
#                 (rowSums(norm_expr[, group_target]) + 1)/length(group_target), 
#                 base = 2)
#             i <- which(names(DEG) == condition_2[1])
#             avg_log2FC <- (mean_target-mean_base)[rownames(DEG[[i]])]
#             output <- data.frame("p_val" = DEG[[i]]$p.value,
#                 "avg_log2FC" = avg_log2FC,
#                 "pct.1" = DEG[[i]]$self.detected,
#                 "pct.2" = DEG[[i]]$other.detected,
#                 "p_val_adj" = DEG[[i]]$FDR)
#             rownames(output) <- rownames(DEG[[i]])
#             output
#         })
#         names(output) <- condition_1
#     }
#     return(output)
# }

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
#' file <- system.file("extdata", "example_mapping_organism.rda",
#' package = "singIST")
#' load(file)
#' data_organism <- example_mapping_organism
#' file <- system.file("extdata", "example_superpathway_fit_model.rda",
#' package = "singIST")
#' load(file)
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
#' \link{diff_expressed} with the logFC for each gene and cell type. 
#' P-values will be corrected with BH, selecting only the genes in the orthologs
#' gene set.
#' @returns
#' A list object with the singIST treated samples predictor block matrix and
#' a list of Fold Changes for each cell type used to compute the singIST
#' treated samples.
#' @import stats
#' @export
#' @examples
#' # Orthology mapping
#' file <- system.file("extdata", "example_mapping_organism.rda",
#' package = "singIST")
#' load(file)
#' data_organism <- example_mapping_organism
#' file <- system.file("extdata", "example_superpathway_fit_model.rda",
#' package = "singIST")
#' load(file)
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
    names(logFC) <- names(object@celltype_mapping)[cells] # gsub("_", " ", names(logFC))
    FC <- vector("list", length(cells))
    for(b in cells){
        genes <- base::intersect(
            orthologs[[b]]$output_gene, rownames(object@counts))
        c <- names(object@celltype_mapping)[b]
        if(length(genes) == 0){
            FC[[c]] <- data.frame()
            next
            }
        FC_aux <- logFC[[c]][rownames(logFC[[c]]) %in% genes, , drop = FALSE]
        FC_aux$`p_val_adj` <- stats::p.adjust(FC_aux[,"p_val"], # Correction pval
                                              method="BH")
        significant_genes <- FC_aux[, "p_val_adj"] <= 0.05
        print(c)
        print(significant_genes)
        if(nrow(FC_aux) == 0){
            FC[[c]] <- data.frame()
            next
            }
        if(sum(significant_genes, na.rm = TRUE) == 0){
            indices_match <- match(rownames(FC_aux), orthologs[[b]]$output_gene)
            FC_aux[!significant_genes, "avg_log2FC"] <-
                rep(0, sum(!significant_genes, na.rm = TRUE))
            rownames(FC_aux) <- paste0(
                c, "*", orthologs[[b]][indices_match, ]$input_gene)
            predictor_block <- FCtoExpression(model_object, b, samples,
                                             predictor_block, FC_aux)
            FC_aux <- FC_aux[, c("avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
            colnames(FC_aux)[1] <- "r_g^b"
            FC[[c]] <- FC_aux
            next
        }
       # FC_aux[significant_genes, "avg_log2FC"] <-
       #        sign(FC_aux[significant_genes, "avg_log2FC"])*
       #    (2^(abs(FC_aux[significant_genes, "avg_log2FC"]))-1)
        FC_aux[!significant_genes, "avg_log2FC"] <-
            rep(0, sum(!significant_genes, na.rm = TRUE))
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
#' @param FC_list Optional parameter with list of matrices containing Fold
#' Changes provided by the user. If such list is provided, logFC will not be 
#' computed via \link{diff_expressed} and the list input will be used instead.
#' The list length must match the number and order of human cell types modelled.
#' Each element of the list must be a `data.frame` whose columns are:
#' "p_value" p-value of test, "avg_log2FC" the log2FC provided, "pct.1" percent
#' of cells where the gene is expressed in base class, "pct.2" percent of cells
#' where the gene is expressed in target class, "p_val_adj" adjusted p-value.
#' Rownames should contain the gene names. Note that log2FC provided should be
#' comparable to log2FC computed in asmbPLS-DA model. log2FC should be reported
#' as descriptive point estimates of mean difference of log2 normalized
#' expression values.
#' @param ... Other parameters to pass onto \link{diff_expressed}
#' @import checkmate SeuratObject
#' @returns
#' A list with; ortholog gene sets as returned by \link{orthology_mapping};
#' a list with the Fold Changes used; singIST treated samples as returned by
#' \link{singIST_treat}
#' @export
biological_link_function <- function(
        object, model_object, object_gene_identifiers = "external_gene_name",
        model_species = "hsapiens", FC_list = NULL, ...){
    # Cell type and orthology mapping
    message("Cell type mapping...")
    object <- celltype_mapping(object)
    #object@counts$test <- paste0(object@counts$celltype_cluster, "_",
    #                               object@counts$class)
    #if(is(object@counts,"Seurat")){SeuratObject::Idents(object@counts)<-"test"}
    if(is.null(FC_list)){
        message("Computing log Fold Changes...")
        logFC <- diff_expressed(object, ...)
    }else{
        checkmate::assert_list(
            FC_list, len = length(names(object@celltype_mapping)))
        checkmate::assert_true(all(names(FC_list) == 
                            gsub("_", " ", names(object@celltype_mapping))))
        for(i in seq_along(FC_list)){
            checkmate::assert_true(all(colnames(FC_list[[i]]) == c("p_val",
                                    "avg_log2FC","pct.1","pct.2","p_val_adj")))
            checkmate::assert_true(all(FC_list[[i]]$p_val_adj <= 1))
        }
        logFC <- FC_list
    }
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
            from_species = model_species)
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