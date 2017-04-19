# File for creating compMS documentation with roxygen

#' Bayesian Compositional Proteomics Modeling via Stan
#'
#' @docType package
#' @name compMS-package
#' @aliases compMS
#' @useDynLib compMS, .registration = TRUE
#'
#' @import methods
#' @importFrom rstan optimizing sampling vb constrain_pars extract extract_sparse_parts get_posterior_mean stanc
#' @import stats
#' @import Rcpp
#' @import bayesplot
#' @import rstantools
#'
#' @description
#' The \pkg{compMS} package enables users to fit compostional proteomics
#' models for isobaric tag mass spectrometry experiments.  The package is
#' designed to make these models accesible to people without a background
#' in Bayesian statistics or compositional data analysis.  This is done by
#' requiring users to follow a very specific format for data entry.  Once
#' the data is formatted properly everything else is handled internally.
#' The required fields and headers are explained in the documentation for the
#' sample data \code{\link{ptmDat}}.  Two sample files are included in the
#' package, the dataframe "sampleDat" and the list object "ptmDat".
#'
#' Once data has been correctly formatted and loaded into your workspace, you can
#' create a Stan model with the function \code{\link{compCall}}.
#'
#' @section Model Details:
#' This is where I will explain the model and the options for use:
#' \describe{
#' Here is a detailed description.
#' }
#'
#'
NULL

