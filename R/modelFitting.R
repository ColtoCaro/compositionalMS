#File for reading in the dataset and fitting the corresponding model

#' Fitting a compositional Stan model
#'
#' This function runs the main code for fitting an appropriate compositonal
#' proteomics model based on the
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
compCall <- function(dat, approx = TRUE, resultsOnly = TRUE, pp=.95){

  #Put single dataframe into a list so that we will always work with a list of dataframes
  if(is.data.frame(dat)){dat <- list(dat)}

  #test to make sure each list component is a dataframe
  testDf <- lapply(dat, is.data.frame)
  if(sum(unlist(testDf)) != length(dat)){
    stop("Error: at least one list component is not a dataframe")
    }

  readyDat <- lapply(1:length(dat), function(x) transformDat(dat[[x]], modelFit,
                                                             x))
  oneDat <- do.call(rbind, readyDat)
  oneDat <- oneDat[order(oneDat$techID, oneDat$ptm, oneDat$ptmID, oneDat$bioID,
                         oneDat$condID, oneDat$tag_plex), ]

  #set data variables
  N_ <- nrow(oneDat)
  n_c <- length(unique(oneDat$techID))

  sumBio <- sum(unlist(lapply(dat, function(x) (x[1, 3] ==1 | x[2,1] == 1)  )))
  if(sumBio == 0){
    n_b <- 0
    bioID <- rep(0,N_)
  }else{
    n_b <- length(unique(oneDat$bioID))
    bioID <- as.integer(factor(oneDat$bioID))
    }

  sumCond <- sum(unlist(lapply(dat, function(x) (x[3, 1] ==1 ))))
  if(sumCond == 0){
    n_gc <- 0
    condID <- rep(0, N_)
  }else{
    n_gc <- length(unique(oneDat$condID))
    condID <- as.integer(factor(oneDat$condID))
    }

  n_t <- length(unique(oneDat$tag_plex))

  sumPtm <- sum(unlist(lapply(dat, function(x) (x[4, 1] == 1))))
  if(sumPtm == 0){
    n_p <- 0
    n_ptm <- 0
    ptm <- rep(0,N_)
    ptmPep <- rep(0,N_)
  }else{
  nonPtms <- which(oneDat$ptm == 0)
  n_p <- length(unique(oneDat[-nonPtms , ]$ptmID))
  n_ptm <- length(unique(oneDat[-nonPtms , ]$ptm))
  ptm <- as.integer(oneDat$ptm)
  ptmPep <- rep(0,N_)
  ptmPep[-nonPtms] <- as.integer(factor(oneDat[-nonPtms , ]$ptmID))
  }

  techID <- as.integer(factor(oneDat$techID))
  tag <- as.integer(factor(oneDat$tag_plex))

  sumCov <- sum(unlist(lapply(dat, function(x) x[1, "Covariate"])))
  useCov <- 1*(sumCov > 0)

  covariate <- oneDat$covariate
  lr <- oneDat$lr

  model <- rstan::sampling(compMS:::stanmodels$allModels)
  model

} #end of compFit function

#' Fitting the model
fitModel <- function(dat, modelFit, approx, resultsOnly){
  #create unique protein id's.
  dat <- lapply(dat, addIds)
}





