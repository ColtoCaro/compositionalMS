#File to test out package functionality

#' @name test_function
#' @param text
#' @export
#' @return Some random text
test_function <- function(text){
  print(paste("You printed:", text))
}

#Can I get at the Stan function?
#' @name callStan
#' @export
cppcode <- rstan::stanc(file = "~/exec/testModel.stan")
