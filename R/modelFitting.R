#File for reading in the dataset and fitting the corresponding model

#' Fitting a compositional Stan model
#'
#' This function fits an appropriate compositonal proteomics model based on the
#' input data file.  It returns a Stan modelfit object.
#'
#' @export
#' @param dat The data to be analyzed.  This data must be formatted as shown in
#'   the included sample data.
#' @param approx A boolean variable which specifies whether or not a full Bayesian
#'   or an approximate Bayesian model should be fit.  The default value of TRUE
#'   fits a Variational Bayes approximation.
#' @details There are many details.  This will be filled out later.
#'
#'
compFit <- function(dat, approx = TRUE){
  NULL
}
