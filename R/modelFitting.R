#File for reading in the dataset and fitting the corresponding model

#' Fitting a compositional Stan model
#'
#' This function fits an appropriate compositonal proteomics model based on the
#' input data file.  It returns a Stan modelfit object.
#'
#' @export
#' @param dat The data to be analyzed.  This data must be formatted as shown in
#'   the included sample data.  Both the reference channel and the reference
#'   conditions will be defined from the first column of data. If multiple runs
#'   are to be analyzed at once then the dat object must be a list of
#'   dataframes.  This is only recommended if biological or technical
#'   replicates exist across multiple runs.  Otherwise, for computational
#'   simplicity, each run should be analyzed separately.  For a single run, the
#'   data should be formatted as a dataframe.  In each dataframe the first
#'   two rows are used to determine which columns are technical replicates and
#'   which are biological replicates.  For example, if two columns, possibly
#'   from different plexes, both have the number '3' in the first row, then
#'   they will be treated as technical replicates.  Likewise, columns that
#'   share a number in the second row are treated as biological replicates. If
#'   there are no biological replicates in the experiment then the second row
#'   should be all zeroes.  If technical replicates are not also labeled as
#'   biological replicates an error will be generated.
#' @param approx A boolean variable which specifies whether or not a full
#'   Bayesian
#'   or an approximate Bayesian model should be fit.  The default value of TRUE
#'   fits a Variational Bayes approximation.
#' @param resultsOnly A boolean variable that determines what information will
#'   be returned after model fitting.  If resultsOnly is TRUE, then the
#'   function will return a list containing two tables, one table for posterior
#'   means and another with posterior standard deviations.  If resultsOnly is
#'   FALSE then the results list will also contain an object of class Stanfit.
#'   The Stanfit object contains all of the sampling and some measures for
#'   quality control.
#' @details There are many details.  This will be filled out later.
#'
#'
compFit <- function(dat, approx = TRUE, resultsOnly = TRUE){
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

  #Check for biological replicates and multiple conditions
  testBioreps <- lapply(dat, function(x) sum(x[2, 4:length(x)],
                                                  na.rm = T))
  if(testBioreps != 0){
    bioReps == TRUE
    }else{bioReps == FALSE}

  testConditions <- lapply(testList, function(x) sum(x[3, 4:length(x)],
                                                  na.rm = T))
  if(testConditions != 0){
    multiCond == TRUE
  }else{multiCond == FALSE}

  #Call the model fitting function determined by multiCond and bioReps
  if(multiCond + bioReps == 0){modelFit <- fitSimple(dat, approx, resultsOnly)}
  if(multiCond + bioReps == 2){modelFit <- fitFull(dat, approx, resultsOnly)}
  if(multiCond + bioReps == 1){
    if(multiCond == TRUE){modelFit <- fitCondition(dat, approx, resultsOnly)}
    if(bioReps == TRUE){modelFit <- fitBioreps(dat, approx, resultsOnly)}
    }

  modelFit

} #end of compFit function

#' Fitting the simple model
#'
#' Function for specifically calling a model with no biological
#' replicates or condition groups.
fitSimple <- function(dat, approx, resultsOnly){
  #create unique protein id's.
  dat <- lapply(dat, addIds)
}




