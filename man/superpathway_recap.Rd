% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step4-recapitulations.R
\name{superpathway_recap}
\alias{superpathway_recap}
\title{Derive superpathway recapitulation}
\usage{
superpathway_recap(model_object, data_original, data_singIST)
}
\arguments{
\item{model_object}{A \link{superpathway.fit.model-class} object passed
from \link{singISTrecapitulations}}

\item{data_original}{A matrix with the superpathway's score as returned
by \link{derive_contributions} for the non-singIST treated samples,
passed from \link{singISTrecapitulations}}

\item{data_singIST}{A matrix with the superpathway's score as returned
by \link{derive_contributions} for the singIST treated samples,
passed from \link{singISTrecapitulations}}
}
\value{
An object \code{data.frame} with the variables: \code{pathway} name as indicated in
\code{model_object}, \code{recapitulation} with the superpathway recapitulation
}
\description{
Derive superpathway recapitulation
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
superpathway_recap(model, original$superpathway_score,
derived$superpathway_score)
}
