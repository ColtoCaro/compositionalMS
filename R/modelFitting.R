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
compCall <- function(dat, approx = FALSE, resultsOnly = FALSE, pp=.95){

  #Put single dataframe into a list so that we will always work with a list of dataframes
  if(is.data.frame(dat)){dat <- list(dat)}

  #test to make sure each list component is a dataframe
  if (length(dat) > 1){
  testDf <- lapply(dat, is.data.frame)
  if(sum(unlist(testDf)) != length(dat)){
    stop("Error: at least one list component is not a dataframe")
  }
  }
  
  #make sure that each dataframe has the same reference condition
  refList <- lapply(dat, function(x) paste(x[1, "tag1"], 
                                           x[2, "tag1"]))
  if(!do.call(all.equal,refList)){
    stop("Error: Plexes have different reference channels")
  }
  
  readyDat <- lapply(1:length(dat), function(x) transformDat(dat[[x]], modelFit,
                                                             x))
  oneDat <- do.call(rbind, readyDat)
  oneDat <- oneDat[order(oneDat$condID, oneDat$bioID, oneDat$ptm, oneDat$ptmID),]
  
  #set data variables
  N_ <- nrow(oneDat)
  n_c <- length(unique(oneDat$condID))
  condKey <- data.frame(number = 1:n_c, name = unique(oneDat$condID))
  condID <- as.integer(factor(oneDat$condID))
  
  sumBio <- sum(unlist(lapply(dat, function(x) (x[1, 3] ==1 | x[2,1] == 1)  )))
  if(sumBio == 0){
    n_b <- 0
    bioID <- rep(0,N_)
    condToBio <- rep(0, n_b)
    n_nc <- rep(0, n_c)
    max_nc <- 0
  }else{
    n_b <- length(unique(oneDat$bioID))
    bioID <- as.integer(factor(oneDat$bioID))
    #make a mapping for use in a heierarchical model (not implemented)
    bioMap <- oneDat$condID[match(unique(oneDat$bioID), oneDat$bioID)]
    conditionNumber <- condKey$number[match(bioMap, condKey$name)]
    bioKey <- data.frame(number = 1:n_b, bioID = unique(oneDat$bioID),
                         condID = bioMap, conditionNumber)
    bioToCond <- bioKey$conditionNumber
    
    #make an array giving bio positions for each condition
    
    n_nc <- unlist(lapply(1:n_c, function(x) sum(bioToCond == x)))
    max_nc <- max(n_nc)
    ncMap <- lapply(1:n_c, function(x) which(bioToCond == x))
    condToBio <- matrix(0, nrow = n_c, ncol = max_nc)
    for (i in 1:n_c){
      condToBio[i, 1:n_nc[i]] <- ncMap[[i]]
    }
  }
  
  
  sumPtm <- sum(unlist(lapply(dat, function(x) (x[3, 1] == 1))))
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
    
    #find and remove ptm data with no corresponding protein data
    globalProts <- unique(oneDat$condID[oneDat$ptm == 0])  
    ptmProts <- unique(oneDat$condID[oneDat$ptm > 0])  
    orphanProts <- setdiff(ptmProts, globalProts)
    if(length(orphanProts > 0)){
      orphanIndex <- which(oneDat$condID %in% orphanProts)
      oneDat <- oneDat[-orphanIndex, ]
      wText <- paste(length(orphanIndex), "PTM data points were removed because they had no corresponding protein level measurements")
      warning(wText)
    }
  } # end actions for ptm experiments
  
  condID <- as.integer(factor(oneDat$condID))
  
  sumCov <- sum(unlist(lapply(dat, function(x) x[1, "Covariate"])))
  useCov <- 1*(sumCov > 0)
  
  covariate <- oneDat$covariate/max(oneDat$covariate)
  lr <- oneDat$lr
  

  sMod <- compMS:::stanmodels$allModels
  if(approx){
    model <- rstan::vb(sMod)
  }else{
    model <- rstan::sampling(sMod)
  }
  model

} #end of compFit function




