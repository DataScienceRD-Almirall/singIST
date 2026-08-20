## inst/scripts/example_orthology_mapping.R
## ---------------------------------------------------------------------------
## Documentation script for inst/extdata/example_orthology_mapping.rda
##
## Describes how the shipped orthology cache was generated. This file is NOT
## run during R CMD check or package build; it only documents provenance so
## that users can reproduce (or refresh) the cached mapping.
##
## Cache contents: 17,066 human -> mouse one-to-one ortholog pairs with
## columns: input_gene, output_gene, input_ensembl, output_ensembl.
##
## It is consumed by orthology_mapping(..., orthology_cache = <this object>)
## so that the vignette and user analyses complete even when Ensembl BioMart
## is unavailable (service interruption or platform migration).
##
## NOTE: record the Ensembl release used when (re)generating, e.g. via
##   biomaRt::listEnsemblArchives()  or  mart@host  after useEnsembl().
## ---------------------------------------------------------------------------

library(biomaRt)

## 1. Connect to Ensembl (human reference + mouse model)
mart_h <- biomaRt::useEnsembl(biomart = "genes",
                              dataset = "hsapiens_gene_ensembl")
mart_m <- biomaRt::useEnsembl(biomart = "genes",
                              dataset = "mmusculus_gene_ensembl")

## 2. All human genes: gene symbol + Ensembl ID
human <- biomaRt::getBM(
    attributes = c("external_gene_name", "ensembl_gene_id"),
    mart = mart_h)
colnames(human) <- c("input_gene", "input_ensembl")

## 3. Mouse one-to-one orthologs for every human Ensembl ID
ortho <- biomaRt::getBM(
    attributes = c("ensembl_gene_id",
                   "mmusculus_homolog_ensembl_gene",
                   "mmusculus_homolog_orthology_type"),
    filters = "ensembl_gene_id",
    values = human$input_ensembl,
    mart = mart_h)
colnames(ortho) <- c("input_ensembl", "output_ensembl", "orthology_type")
ortho <- ortho[ortho$orthology_type == "ortholog_one2one", ]

## 4. Mouse Ensembl ID -> gene symbol
mouse <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    mart = mart_m)
colnames(mouse) <- c("output_ensembl", "output_gene")

## 5. Assemble the cache in the expected column order
example_orthology_mapping <- merge(human, ortho, by = "input_ensembl")
example_orthology_mapping <- merge(example_orthology_mapping, mouse,
                                   by = "output_ensembl")
example_orthology_mapping <- example_orthology_mapping[
    , c("input_gene", "output_gene", "input_ensembl", "output_ensembl")]

## 6. Save (bzip2 compression keeps the file well under the 5 MB per-file
##    Bioconductor limit; the shipped file is ~168 KB)
save(example_orthology_mapping,
     file = "inst/extdata/example_orthology_mapping.rda",
     compress = "bzip2")
