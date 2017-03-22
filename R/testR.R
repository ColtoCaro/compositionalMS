#File to test out package functionality

#' Testing this stuff out
#'
#' This code is just a way for me to learn roxygen and to test out
#' access to Stqn objects
#'
#' @export
#' @param y_ an input vector
#'
#' @return The posterior mean of beta
#'
#' @details We fit the model y = beta + epsilon.  Hopefully we will return
#'   the correct mean.
#' @example
#' test_function(rnorm(100, 5, 3))
test_function <- function(y_){
  N <- length(y_)
  tempMod <- sampling(stanmodels$testModel)
  tempMod
}


