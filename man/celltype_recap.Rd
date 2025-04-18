% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step4-recapitulations.R
\name{celltype_recap}
\alias{celltype_recap}
\title{Derive cell type recapitulation}
\usage{
celltype_recap(model_object, data_original, data_singIST)
}
\arguments{
\item{model_object}{A \link{superpathway.fit.model-class} object passed
from \link{singISTrecapitulations}}

\item{data_original}{A matrix with the cell type contributions as returned
by \link{derive_contributions} for the non-singIST treated samples,
passed from \link{singISTrecapitulations}}

\item{data_singIST}{A matrix with the cell type contributions as returned
by \link{derive_contributions} for the singIST treated samples,
passed from \link{singISTrecapitulations}}
}
\value{
A \code{data.frame} object with the variables: \code{pathway} name, \code{celltype} with
the cell type name, \code{recapitulation} with the cell type recapitulation, and
\code{reference} with the cell type reference recapitulation
}
\description{
Derive cell type recapitulation
}
\examples{
data(example_superpathway_fit_model)
model <- example_superpathway_fit_model
data(example_mapping_organism)
mapped <- example_mapping_organism
singIST_samples <- biological_link_function(mapped, model)$singIST_samples
original <- derive_contributions(model, singIST_samples)
derived <- derive_contributions(model,
slot(model, "model_fit")$predictor_block)
celltype_recap(model, original$celltype_contribution,
derived$celltype_contribution)
}
