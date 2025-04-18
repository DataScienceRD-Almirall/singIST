% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbb-step1-fit.R
\name{asmbPLSDA.cv.loo}
\alias{asmbPLSDA.cv.loo}
\title{Leave-one-out Cross-validation}
\usage{
asmbPLSDA.cv.loo(
  X.matrix,
  Y.matrix,
  PLS_term = 1,
  X.dim,
  quantile.comb.table,
  outcome.type = c("binary", "multiclass"),
  Method = NULL,
  measure = "B_accuracy",
  parallel = FALSE,
  expected.measure.increase = 0.005,
  center = TRUE,
  scale = TRUE,
  maxiter = 100
)
}
\arguments{
\item{X.matrix}{Predictor block matrix from \code{matrixToBlock}}

\item{Y.matrix}{Response matrix from \code{matrixToBlock}}

\item{PLS_term}{An integer with the number of PLS components to use passed
from \link{hyperparameters-class} obect}

\item{X.dim}{A list with the observed gene set size for each cell type
from \code{matrixToBlock}}

\item{quantile.comb.table}{A matrix with the quantile comb table passed
from \link{hyperparameters-class} object}

\item{outcome.type}{A character indicating \code{binary} or \code{multiclass} passed
from \link{hyperparameters-class} object}

\item{Method}{A parameter passed from \code{fitOptimal}}

\item{measure}{A parameter passed from \code{fitOptimal}}

\item{parallel}{A parameter passed from \code{fitOptimal}}

\item{expected.measure.increase}{A parameter passed from \code{fitOptimal}}

\item{center}{A parameter passed from \code{fitOptimal}}

\item{scale}{A parameter passed from \code{fitOptimal}}

\item{maxiter}{A parameter passed from \code{fitOptimal}}
}
\value{
A list containing the optimal quantiles for each PLS component and the
optimal number of PLS components.
}
\description{
Leave-one-out Cross-validation
}
\examples{
data(example_superpathway_input)
data <- example_superpathway_input
matrices <- matrixToBlock(data)
X.matrix <- matrices$block_predictor
Y.matrix <- matrices$matrix_response
X.dim <- matrices$block_dim
quantile.comb.table <- slot(slot(data, "hyperparameters_info"),
"quantile_comb_table")
outcome.type <- slot(slot(data, "hyperparameters_info"), "outcome_type")
asmbPLSDA.cv.loo(X.matrix, Y.matrix, PLS_term = 1, X.dim,quantile.comb.table,
Method = NULL, measure = "B_accuracy", parallel = TRUE,
outcome.type = outcome.type, expected.measure.increase = 0.005,
center = TRUE, scale = TRUE,maxiter = 100)
}
