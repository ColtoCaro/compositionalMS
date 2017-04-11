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
compCall <- function(dat, approx = TRUE, resultsOnly = TRUE){

  #Put single dataframe into a list so that we will always work with a list of dataframes
  if(is.data.frame(dat)){dat <- list(dat)}

  #test to make sure each list component is a dataframe
  testDf <- lapply(dat, is.data.frame)
  if(sum(unlist(testDf)) != length(dat)){
    stop("Error: at least one list component is not a dataframe")
    }

  #make sure column names are valid
  if(checkQs(dat) > 0){stop("Column names cannot contain qqqq")}
  if(checkNames(dat) > 0){stop("Column names must be unique")}

  #Check for biological replicates, multiple conditions and covariates
  sumBioreps <- sum(unlist(lapply(dat, function(x) sum(x[2, 4:length(x)],
                                                  na.rm = T))))
  if(sumBioreps != 0){
    bioReps <- TRUE
    }else{bioReps <- FALSE}

  sumConditions <- sum(unlist(lapply(dat, function(x) sum(x[3, 4:length(x)],
                                                  na.rm = T))))
  if(sumConditions != 0){
    multiCond <- TRUE
  }else{multiCond <- FALSE}

  testCovariate <- sum(unlist(lapply(dat, function(x) x[1, 3])))
  if(testCovariate == 0){covariate <- FALSE}
  if(testCovariate == length(dat)){covariate <- TRUE}
  if(!(testCovariate == 0 | testCovariate == length(dat))){
    stop("error: inconsistent covariate use")}

  #Call the model fitting function determined by multiCond, and bioReps
  if(multiCond + bioReps == 0){
    if(covariate == TRUE){modelFit <- "simpleSlope"}
    if(covariate == FALSE){modelFit <- "simpleInt"}
  }
  if(multiCond + bioReps == 2){
    if(covariate == TRUE){modelFit <- "fullSlope"}
    if(covariate == FALSE){modelFit <- "fullInt"}
  }
  if(multiCond + bioReps == 1){
    if(multiCond == TRUE){
      if(covariate == TRUE){modelFit <- "condSlope"}
      if(covariate == FALSE){modelFit <- "condInt"}
    }
    if(bioReps == TRUE){
      if(covariate == TRUE){modelFit <- "biorepSlope"}
      if(covariate == FALSE){modelFit <- "biorepInt"}
    }
    }

  readyDat <- transformDat(dat, modelFit)
  model <- fitModel(readyDat, modelFit, approx, resultsOnly)
  model

} #end of compFit function

#' Fitting the model
fitModel <- function(dat, modelFit, approx, resultsOnly){
  #create unique protein id's.
  dat <- lapply(dat, addIds)
}


gProt <- unique(global$ProteinId)
phosInt <- phos[phos$ProteinId %in% gProt, ]
ubiqInt <- ubiq[ubiq$ProteinId %in% gProt, ]

globalSSN <- rowSums(global[ , grep("rq", colnames(global))])
phosSSN <- rowSums(phosInt[ , grep("rq", colnames(phosInt))])
ubiqSSN <- rowSums(ubiqInt[ , grep("rq", colnames(ubiqInt))])

protDat <- data.frame(Protein = global$ProteinId, Peptide = global$PeptideSequence,
                      bioID = global$Class, Covariate = globalSSN,
                      global[ , grep("rq", colnames(global))])
phosDat <- data.frame(Protein = phosInt$ProteinId, Peptide = phosInt$PeptideSequence,
                      bioID = phosInt$Class, Covariate = phosSSN,
                      phosInt[ , grep("rq", colnames(phosInt))])
ubiqDat <- data.frame(Protein = ubiqInt$ProteinId, Peptide = ubiqInt$PeptideSequence,
                      bioID = ubiqInt$Class, Covariate = ubiqSSN,
                      ubiqInt[ , grep("rq", colnames(ubiqInt))])

protHead <- matrix(0, ncol = 14, nrow= 4)
rownames(protHead) <- c("techReps", "bioReps", "Condition", "PTM")
protHead[1, ] <- c(1, 0,1, 1, 1:10)

phosHead <- protHead
ubiqHead <- protHead

phosHead[4, ] <- c(1, 0, 0, 0, rep(1,10))
ubiqHead[4, ] <- c(1, 0, 0, 0, rep(2,10))

protHead <- as.data.frame(protHead)
phosHead <- as.data.frame(phosHead)
ubiqHead <- as.data.frame(ubiqHead)

names(protHead) <- names(phosHead) <- names(ubiqHead) <- names(protDat)
ptmDat1 <- rbind(protHead, protDat)
ptmDat2 <- rbind(phosHead, phosDat)
ptmDat3 <- rbind(ubiqHead, ubiqDat)

save(ptmDat1, file="/Users/darkapple/Documents/compMS/data/ptmDat1.Rdata")
save(ptmDat2, file="/Users/darkapple/Documents/compMS/data/ptmDat2.Rdata")
save(ptmDat3, file="/Users/darkapple/Documents/compMS/data/ptmDat3.Rdata")

