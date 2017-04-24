#File for functions that create plots

#' Caterpillar plot
#'
#' This function takes a list which contains a Stanfit object and summary
#' tables.  The output is a caterpillar plot of parameters with credible
#' intervals that do not contain 0.  Red line segments denote a 80% credible
#' interval, while the black tails show 95%.
#'
#' @export
#' @param results The results list generated from the function
#'   \code{\link{compCall}}
#' @param ptm A boolean variable denoting whether or not plots for ptm peptides
#'   should be generated.
#' @param allPars A boolean variable indicating whether or not all parameters
#'   should be plotted.  The default is false, resulting in a plot for only
#'   parameters with 95% credible intervals that do not contain zero.
#'
caterpillar <- function(results, ptm = FALSE, allPars = FALSE){
  if(ptm == T){
    parSumm <- rstan::summary(results[[3]], pars = "alpha")$summary
    parName <- "alpha"
  }else{
    if(results[[3]]@par_dims$avgCond == 0){
      parSumm <- rstan::summary(results[[3]], pars = "beta")$summary
      parName <- "beta"
    }else{
      parSumm <- rstan::summary(results[[3]], pars = "avgCond")$summary
      parName <- "avgCond"
    }
  }

  if(allPars == FALSE){
    sigParIndex <- which(parSumm[ , "2.5%"] > 0 | parSumm[ , "97.5%"] < 0)
  }else{
    sigParIndex <- 1:dim(parSumm)[1]
  }
  orderedIndex <- sigParIndex[order(parSumm[sigParIndex, "50%"])]

  parStr <- paste(parName, "[", orderedIndex, "]", sep = "")

  rstan::plot(results[[3]], pars = parStr, mapping = theme) +
    scale_y_discrete(labels = NULL)

}#end caterpillar

#' Precision plot
#'
#' This function takes a results summary from the \code{\link{compCall}}
#' function.
#' The output is a plot of showing posterior means on the
#' x-axis and precision on the y-axis, where precision defined as the inverse
#' of the posterior variance.
#'
#' @export
#' @param summary A results table generated from the function
#'   \code{\link{compCall}}, which will be found in either the first (for
#'   protein summaries) or second (for ptm's) component of the result list
#'   object.
#'
precisionPlot <- function(summary){
  ggplot2::ggplot(results[[1]], aes(x = mean, y = 1/var)) +
    ggplot2::geom_point(aes(color = P_null)) + scale_color_gradient(low = "red", high = "black")
}

#'
#' Variance plot
#'
#' This function shows the distribution of each parameter corresponding to an
#' experimental error.  Different variance components are modeled for each
#' tag/plex combination (along with different components for each PTM). If one
#' parameter is greatly different from the others it could be indicative of
#' problems in the experimental workflow.
#'
#' @export
#' @param model A Stanmodel object generated from the function
#' \code{\link{compCall}}, which can be found in the third component of the
#' results object.
#'
checkVariance <- function(model){
  rstan::plot(model, pars ="sigma")
}

