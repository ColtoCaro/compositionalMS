 #File for reading in the dataset and fitting the corresponding model

#' Fitting a compositional Stan model
#'
#' This function runs the main code for fitting an appropriate compositonal
#' proteomics model based on the
#' input data file.  It returns a list of size 3.  The first component of the
#' contains summary data on relative protein abundance.  The second component
#' summarizes PTM peptides and the third component contains the Stanfit object.
#'
#' @export
#' @param dat The data to be analyzed.  This data must be formatted as shown in
#'   the included sample data \code{\link{ptmDat}}.  Both the reference channel
#'    and the reference
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
#' @param pp The percentage covariate level to be used for predicting relative
#'   protein abundance.  By default this value is .95, so that sum signal to
#'   noise will be adjusted towards the 95th percentile of observed values.
#' @param nCores determines the number of cores used for parallel processing
#'   should be used.
#'   For large datasets it is highly recommended that the package be used on a
#'   computer capable of parallel processing.
#' @param iter Number of iterations for each chain
#' @param normalize A boolean variable that determines whether or not
#'   adjustments should be made under the assumption that the average protein
#'   abundance is equivalent in each channel.  By default this value is true
#'   which results in each row of the matrix being perturbed by the inverse of
#'   the compositional mean.  Consequently the geometric means of each tag will
#'   be equivalent.
#'
#' @details There are many details.  This will be filled out later.
#'
#'
#'
compBayes <- function(dat,
                     pp=.99,
                     nCores = 1,
                     iter = 2000,
                     normalize = TRUE
                     ){

  #Put single dataframe into a list so that we will always work with a list of dataframes
  if(is.data.frame(dat)){dat <- list(dat)}

  #test to make sure each list component is a dataframe
  if (length(dat) > 1){
    testDf <- lapply(dat, is.data.frame)
    if(sum(unlist(testDf)) != length(dat)){
      stop("Error: at least one list component is not a dataframe")
    }
    #make sure that each dataframe has the same reference condition
    refList <- lapply(dat, function(x) paste(x[1, "tag1"],
                                             x[2, "tag1"]))
    if(!do.call(all.equal,refList)){
      stop("Error: Plexes have different reference channels")
    }
  }

  #make sure protein names do not have "_" characters
  dat <- lapply(dat, function(x) within(x, Protein <-
                  gsub("_", "-", Protein)))

  #determine how many redundancies are being used
  maxRedun <- max(unlist(lapply(dat, function(x) max(x[1, "varCat"]))))


  readyDat <- lapply(1:length(dat), function(x)
    transformDat(dat[[x]], plexNumber = x, normalize = normalize))
  oneDat <- do.call(rbind, readyDat)

  #set data variables
  #Do ptms first since it might change the dimensions of the data
  N_ <- nrow(oneDat)
  sumPtm <- sum(unlist(lapply(dat, function(x) (x[3, 1] == 1))))
  if(sumPtm == 0){
    n_p <- 0
    n_ptm <- 0
    ptm <- rep(0,N_)
    ptmPep <- rep(0,N_)
  }else{
    stop("You have requested a PTM analysis.  This feature is currently
         under development.")
    #find and remove ptm data with no corresponding protein data
    globalProts <- unique(oneDat$bioID[oneDat$ptm == 0])
    ptmProts <- unique(oneDat$bioID[oneDat$ptm > 0])
    orphanProts <- setdiff(ptmProts, globalProts)
    if(length(orphanProts > 0)){
      #remove orphan prots from ptm data
      ptmIndex <- which(oneDat$ptm > 0)
      orphanIndex <- which(oneDat$bioID %in% orphanProts)
      bad <- intersect(ptmIndex, orphanIndex)
      oneDat <- oneDat[-bad, ]
      wText <- paste(length(orphanIndex), "PTM data points were removed because they had no corresponding protein level measurements")
      warning(wText)
    }

    nonPtms <- which(oneDat$ptm == 0)
    ptmName <- levels(factor(oneDat[-nonPtms , ]$ptmID))
    n_p <- length(ptmName)
    n_ptm <- length(unique(oneDat[-nonPtms , ]$ptm))
    ptm <- as.integer(oneDat$ptm)
    ptmPep <- rep(0, nrow(oneDat))
    ptmPep[-nonPtms] <- as.integer(factor(oneDat[-nonPtms , ]$ptmID))

  } # end actions for ptm experiments

  oneDat <- oneDat[order(oneDat$condID, oneDat$bioID, oneDat$ptm,
                         oneDat$ptmID),]

  N_ <- nrow(oneDat)
  n_c <- length(unique(oneDat$condID))
  condKey <- data.frame(number = 1:n_c,
                        name = levels(factor(oneDat$condID)))
  condID <- as.integer(factor(oneDat$condID))

  #Create a tag by varCat variable which determines vc's
  if(maxRedun == 0){
    tagID <- as.integer(factor(oneDat$tag_plex))
  }else{
    redunStr <- paste("R", oneDat$varCat, sep = "")
    if(sumPtm > 0){
      redunStr[which(oneDat$ptm > 0)] <- "Rp"
    }
    oneDat$tag_plex <- paste(oneDat$tag_plex, redunStr, sep="")
    tagID <- as.integer(factor(oneDat$tag_plex))
    }


  n_t <- length(unique(oneDat$tag_plex))


  sumBio <- sum(unlist(lapply(dat, function(x) (x[1, 3] ==1 | x[2,1] == 1)  )))
  if(sumBio == 0){
    n_b <- 0
    bioID <- rep(0,N_)
    condToBio <- matrix(0, nrow = n_c, ncol = 1)
    n_nc <- rep(0, n_c)
    max_nc <- 0
  }else{
    stop("You have requested a population level.  This feature is
         currently under development")
    n_b <- length(unique(oneDat$bioID))
    bioID <- as.integer(factor(oneDat$bioID))
    #make a mapping for use in a heierarchical model (not yet implemented)
    bioMap <- oneDat$condID[match(levels(factor(oneDat$bioID)),
                                  oneDat$bioID)]
    conditionNumber <- condKey$number[match(bioMap, condKey$name)]
    bioKey <- data.frame(number = 1:n_b, bioID = levels(factor(oneDat$bioID)),
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


  sumCov <- sum(unlist(lapply(dat, function(x) x[1, "Covariate"])))
  useCov <- 1*(sumCov > 0)

  if(useCov){
    covariate <- oneDat$covariate/quantile(oneDat$covariate, probs = pp)
  }else{covariate <- oneDat$covariate}

  lr <- oneDat$lr
  if(sum(lr) == 0){stop("Outcomes are all zero. This might be the
                        consequence of normalizing values already less than
                        one")}

  summaryStr <- paste("Estimating ", max(n_b, n_c), " relative protein abundances, and ", n_p, "protein adjusted ptm changes")
  print(summaryStr)

  #local call for testing
 # model <- stan(file="~/Documents/compMS/exec/allModels.stan",
   #             iter = 2000, cores = 4, control = list(adapt_delta = .8))

  sMod <- compMS:::stanmodels$allModels
  model <- rstan::sampling(sMod, cores = nCores, iter = iter)


  #create summary table
  if(n_b > 0){
    targetChain <- rstan::extract(model, pars="avgCond")$avgCond
    postMeans <- colMeans(targetChain)
    postVar <- apply(targetChain, 2, var)
    pvals <- pnorm(nullSet[2], postMeans, sqrt(postVar)) -
      pnorm(nullSet[1], postMeans, sqrt(postVar))
    pValue <- 1 - pchisq((postMeans^2)/postVar, 1)
    justProts <- oneDat[oneDat$ptm == 0 , ]
    n_obs <- unlist(by(justProts$lr, justProts$condID, length))
  }else{
    if(useCov == 0){
      targetChain <- rstan::extract(model, pars="beta")$beta
    }else{
      targetChain <- rstan::extract(model, pars="betaP_c")$betaP_c
    }
    postMeans <- colMeans(targetChain)
    postVar <- apply(targetChain, 2, var)
    intervals <- t(apply(targetChain, 2, quantile,
                         probs = c(.025, .975, .1, .9)))
    colnames(intervals) <- c("LL95", "UL95", "LL80", "UL80")
    justProts <- oneDat[oneDat$ptm == 0, ]
    n_obs <- unlist(by(justProts$lr, justProts$condID, length))
  }

  resDf <- data.frame(name = levels(factor(oneDat$condID)), mean = postMeans,
                      var = postVar, intervals, n_obs = as.vector(n_obs),
                      stringsAsFactors = F)

  if(n_p > 0){
    targetChain <- rstan::extract(model, pars="alpha")$alpha
    postMeans <- colMeans(targetChain)
    postVar <- apply(targetChain, 2, var)
    pvals <- pnorm(nullSet[2], postMeans, sqrt(postVar)) -
      pnorm(nullSet[1], postMeans, sqrt(postVar))
    pValue <- 1 - pchisq((postMeans^2)/postVar, 1)
    justPtms <- oneDat[oneDat$ptm > 0 , ]
    n_obs <- unlist(by(justPtms$lr, justPtms$ptmID, length))
    ptmDf <- data.frame(ptmName, mean = postMeans, var = postVar,
                        P_null = pvals, pValue, n_obs = as.vector(n_obs),
                        stringsAsFactors = F)
  }else{
    ptmDf <- NULL
  }


  RES <- list()
  RES[[1]] <- resDf
  RES[[2]] <- ptmDf
  RES[[3]] <- model


  RES
} #end of compBayes function




#' A transformation from log ratio space into proportions
#'
#' This function transforms all estimated log ratios into proportions
#'
#' @export
#' @param RES Results from the compBayes function.
#'
#' @details Estimates from compBayes are done after using the additive log ratio
#'   transformation.  This function applies the inverse transformation to the
#'   full sampling chains and recomputes summary statistics.  Returns a data
#'   frame with estimates, variance and intervals in terms of proportions.
#'
#'
#'
toSimplex <- function(RES){

  condNames <- RES[[1]]$name
  condNum <- getCond(condNames)
  nCond <- length(unique(condNum))
  nProts <- length(condNames)/nCond

  indices <- lapply(1:nCond, function(x) which(condNum == x + 1))
  indices <- do.call(cbind, indices)

  simpRES <- lapply(1:nProts, function(x)
    alrInv(makeMat(indices[x, ], RES[[3]])))

  #now manipulate results into a dataframe
  pName <- substr(condNames[indices[ , 1]], 1,
                  regexpr("_", condNames[indices[ , 1]])-1)
  Condition <- rep(1:(nCond + 1), nProts)
  Name <- rep(pName, each = nCond + 1)

  simpMat <- do.call(rbind, simpRES)
  simpDf <- data.frame(Name, Condition, simpMat)

  simpDf
} #end to Simplex



