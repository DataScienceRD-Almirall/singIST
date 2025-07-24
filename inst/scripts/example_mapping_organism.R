## code to prepare `example_mapping_organism` object goes here
counts <- SeuratObject::LayerData(SeuratObject::pbmc_small, "counts")
# Change gene symbols to mouse gene symbols = lower case human gene symbol
CapStr <- function(y) {
    c <- strsplit(y, " ")[[1]]
    paste(toupper(substring(c, 1,1)), substring(c, 2),
          sep="", collapse=" ")
}
capitalize_str <- function(charcter_string){
    sapply(charcter_string, CapStr)
}
rownames(counts) <- capitalize_str(tolower(rownames(counts)))
# Update meta.data columns to include class and celltype_cluster variables
# as required by "mapping.organism" class
data <- Seurat::CreateSeuratObject(counts)
data@meta.data$class<- SeuratObject::pbmc_small$groups
data@meta.data$celltype_cluster <- SeuratObject::pbmc_small$RNA_snn_res.1
data@meta.data$donor <- SeuratObject::pbmc_small$letter.idents
# Lognormalize counts and create data layer
data <- Seurat::NormalizeData(data, normalization.method = "LogNormalize")
# Initialize mapping.organism object
example_mapping_organism <- new("mapping.organism",
                                organism = "Mus musculus",
                                target_class = "g1",
                                base_class = "g2",
                                celltype_mapping = list("Dendritic Cells" = c(0,2),
                                                        "Keratinocytes" = c(1)),
                                counts = data)
usethis::use_data(example_mapping_organism, overwrite = TRUE)
