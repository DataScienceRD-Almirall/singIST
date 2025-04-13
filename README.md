
<!-- README.md is generated from README.Rmd. Please edit that file -->

# singIST

<!-- badges: start -->
<!-- badges: end -->


![gabst_01dic](https://github.com/user-attachments/assets/0d36443d-007f-423b-ae75-64bb1a3e23c2)

This repository includes an R library to perform an end-to-end singIST analysis for multiple superpathways and disease models of interest. singIST is a method for comparative single-cell transcriptomics between disease models and humans, which provides explainable and quantitative insights into single-cell transcriptomics alignment. 

## Table of Contents

1. [Library structure](#Library-structure)
2. [Installation](#Installation)
3. [Vignette](#Vignette)

# Library structure

The R library code, R folder, is structured in the following manner:
1. aaa-all-classes.R: contains all classes defined in singIST.
2. ccc-all-methods.R: contains all methods defined in singIST.
3. bbb-step1-fit.R: contains all functions to perform step 1 from singIST, which includes: fit and cross validation of optimal model asmbPLS-DA, and validation metrics of the optimal model. The main function of this code is fitOptimal(), which performs all the aforementioned processes.
4. step2-biologicalLink.R: contains all functions to perform step 2 from singIST, which includes: cell type mapping, orthology mapping and translating the fold changes of the disease models to the human single-cell expression. The main function of this code is biological_link_function() this performs all the aforementioned processes. 
5. step3-contributions.R: it contains the functions to derive the superpathway score, cell type contributions to superpathway score, and gene contributions to superpathway score. These functions are necessary to compute the reference and predicted recapitulations. The main function of this code is derive_contributions().
6. step4-recapitulations.R: it contains the functions to compute the recapitulations at superpathway and cell type levels, including also the computation of the gene contributions to the cell type recapitulation. The main function is singISTrecapitulations().
7. zzz-wrappers.R: defines two wrappers, for the functions that are mainly going to be used for the user, one is multiple_fitOptimal() which performs fitOptimal() for multiple superpathways of interest, and multiple_singISTrecapitulations() which computes recapitulations for a disease model across multiple superpathways of interest.
8. helpers.R: secondary and auxiliary functions used in all the former.

The main functions to be used as a user are fitOptimal(), or its wrapper multiple_fitOptimal(), and singISTrecapitulations(), or its wrapper multiple_singISTrecapitulations(). 

# Installation

``` r
remotes::install_github("DataScienceRD-Almirall/singIST")
```

# Vignette

Find the vignette here https://github.com/DataScienceRD-Almirall/singIST/tree/master/vignettes