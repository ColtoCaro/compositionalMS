#File for reading in the dataset and fitting the corresponding model

#' Fitting a compositional Stan model
#'
#' This function fits an appropriate compositonal proteomics model based on the
#' input data file.  It returns a Stan modelfit object.
#'
#' @export
#' @param dat The data to be analyzed.  This data must be formatted as shown in
#'   the included sample data.  If multiple plexes are to be analyzed at once then
#'   the dat object must be a list of dataframes.  For a single run, the data
#'   should be formatted as a dataframe.  In each dataframe the first
#'   two rows are used to determine which columns are technical replicates and
#'   which are biological replicates.  For example, if two columns, possibly
#'   from different plexes, both have the number '3' in the first row, then
#'   they will be treated as technical replicates.  Likewise, columns that
#'   share a number in the second row are treated as biological replicates.  If
#'   technical replicates are not also labeled as biological replicates an
#'   error will be generated.
#' @param approx A boolean variable which specifies whether or not a full
#'   Bayesian
#'   or an approximate Bayesian model should be fit.  The default value of TRUE
#'   fits a Variational Bayes approximation.
#' @details There are many details.  This will be filled out later.
#'
#'
compFit <- function(dat, approx = TRUE){
  #Determine if there are multiple runs. Make sure data is formatted correctly
  if(is.data.frame(dat) == TRUE){multiRun <- False}
  if(is.list == TRUE){
    #test to make sure each list component is a dataframe
    testDf <- lapply(dat, is.data.frame)
    if(sum(unlist(testDf)) != length(dat)){
      print("Error: at least one list component is not a dataframe")
      return("Error")
    }
  }else{multiRun <- TRUE}

  #Put single dataframe into a list so that we will always work with a list of dataframes
  if(multiRun == FALSE){dat <- list(dat)}

  #Check for biological replicates



} #end of compFit function

testList <- list(sampleDat, sampleDat)
