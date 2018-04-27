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
#'   the included sample data \code{\link{sampleDat}}.  Both the reference channel
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
#' @param pop_sd The prior standard deviation of the population level variance.
#' @param simpleMod A boolean variable that determines whether or not a population
#'   level model will be fit.  In a simple model only relative abuandace parameters
#'   that define average log-ratios within sample groups will be estimated.  The
#'   population level model estimates parameters for each biological replicate along
#'   with population level effects.  Average sample parameters are also estimated
#'   as many experiments do not have sufficient sample size to make use of the full
#'   population level model.
#'
#' @details There are many details.  This will be filled out later.
#'
#'
#'
compBayes <- function(dat,
                     pp=.99,
                     nCores = 1,
                     iter = 2000,
                     normalize = TRUE,
                     pop_sd = 10,
                     simpleMod = TRUE
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
    if(!do.call(all.equal, refList)){
      stop("Error: Plexes have different reference channels")
    }
  }


  #make sure protein names do not have "_" characters
  dat <- lapply(dat, function(x) within(x, Protein <-
                  gsub("_", "-", Protein)))

  #determine how many redundancies are being used
  maxRedun <- max(unlist(lapply(dat, function(x) max(x[1, "varCat"]))))


  readyDat <- lapply(1:length(dat), function(x)
    transformDat(dat[[x]], plexNumber = x, normalize = normalize,
                 simpleMod))
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
                         oneDat$ptmID), ]

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


  #sumBio <- sum(unlist(lapply(dat, function(x) (x[1, 3] ==1 | x[2,1] == 1)  )))
  sumBio <- 1
  if(sumBio == 0){
    n_b <- 0
    bioID <- rep(0,N_)
    condToBio <- matrix(0, nrow = n_c, ncol = 1)
    n_nc <- rep(0, n_c)
    max_nc <- 0
  }else{
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

  if(simpleMod){
    summaryStr <- paste("Estimating ", n_c, " relative protein abundances, and ", n_p, "protein adjusted ptm changes")
  }else{
    summaryStr <- paste("Estimating ", max(n_b, n_c), " relative protein abundances, and ", n_p, "protein adjusted ptm changes")
  }
    print(summaryStr)

  #local call for testing
 # model <- stan(file="~/Documents/compMS/exec/allModels.stan",
   #             iter = 2000, cores = 4, control = list(adapt_delta = .8))

  sMod <- compMS:::stanmodels$allModels
  model <- rstan::sampling(sMod, cores = nCores, iter = iter)


  #create summary tables

  #simple case
  if((n_b == n_c) | simpleMod){

      uBio <- levels(factor(oneDat$condID))
      condNum <- getCond(uBio, bio = FALSE)
      condNames <- getName(uBio)

      refC <- dat[[1]][1, "tag1"]

      nCond <- length(unique(condNum))
      uCond <- unique(c(refC, unique(condNum)))
      uCond <- uCond[order(uCond)]

      nProts <- length(unique(condNames))

      indices <- lapply(1:nProts, function(x)
        which(condNames == levels(factor(condNames))[x]))

      condList <- lapply(indices, function(x) condNum[x])

      simpRES <- lapply(indices, function(x)
        alrInv(makeMat(x, model, bio = FALSE, useCov, avgCond = FALSE)))


      refPos <- which(uCond == refC)
      suCond <- uCond[-refPos]
      if(length(suCond) == 1){
        lrs <- as.matrix(mapply(function(x, y) x[[2]]$Estimate[match(suCond, y)], simpRES, condList))
        colnames(lrs) <- paste0("Est_avg_fc", suCond)
        lrVars <- as.matrix(mapply(function(x, y) x[[2]]$Variance[match(suCond, y)], simpRES, condList))
        colnames(lrVars) <- paste0("Var_Cond", suCond)
        lrLL95 <- as.matrix(mapply(function(x, y) x[[2]]$LL95[match(suCond, y)], simpRES, condList))
        colnames(lrLL95) <- paste0("LL95_Cond", suCond)
        lrUL95 <- as.matrix(mapply(function(x, y) x[[2]]$UL95[match(suCond, y)], simpRES, condList))
        colnames(lrUL95) <- paste0("UL95_Cond", suCond)
      }else{
        lrs <- t(mapply(function(x, y) x[[2]]$Estimate[match(suCond, y)], simpRES, condList))
        colnames(lrs) <- paste0("Est_avg_fc", suCond)
        lrVars <- t(mapply(function(x, y) x[[2]]$Variance[match(suCond, y)], simpRES, condList))
        colnames(lrVars) <- paste0("Var_Cond", suCond)
        lrLL95 <- t(mapply(function(x, y) x[[2]]$LL95[match(suCond, y)], simpRES, condList))
        colnames(lrLL95) <- paste0("LL95_Cond", suCond)
        lrUL95 <- t(mapply(function(x, y) x[[2]]$UL95[match(suCond, y)], simpRES, condList))
        colnames(lrUL95) <- paste0("UL95_Cond", suCond)
      }
      avgLrTab <- data.frame(Protein = levels(factor(condNames)), lrs, lrVars, lrLL95, lrUL95,
                             stringsAsFactors = F)

      #Now get the proportions
      lrs <- t(mapply(function(x, y) x[[1]]$Estimate[match(uCond, unique(c(refC, y)))], simpRES, condList))
      colnames(lrs) <- paste0("Est_avg_Prop", uCond)
      lrVars <- t(mapply(function(x, y) x[[1]]$Variance[match(uCond, unique(c(refC, y)))], simpRES, condList))
      colnames(lrVars) <- paste0("Var_Cond", uCond)
      lrLL95 <- t(mapply(function(x, y) x[[1]]$LL95[match(uCond, unique(c(refC, y)))], simpRES, condList))
      colnames(lrLL95) <- paste0("LL95_Cond", uCond)
      lrUL95 <- t(mapply(function(x, y) x[[1]]$UL95[match(uCond, unique(c(refC, y)))], simpRES, condList))
      colnames(lrUL95) <- paste0("UL95_Cond", uCond)

      avgPTab <- data.frame(Protein = levels(factor(condNames)), lrs, lrVars, lrLL95, lrUL95,
                            stringsAsFactors = F)


      bioDf <- NULL
      popLrTab <- NULL
      popPTab  <- NULL

  }else{    #end simple case

  #Make Biorep proportion table
  uBio <- levels(factor(oneDat$bioID))
  condNum <- getCond(uBio, bio = TRUE)
  condNames <- getName(uBio)

  refC <- dat[[1]][2, "tag1"]
  nCond <- length(unique(condNum)) + 1
  uCond <- unique(c(refC, unique(condNum)))
  uCond <- uCond[order(uCond)]

  nProts <- length(unique(condNames))
  indices <- lapply(1:nProts, function(x) which(condNames == levels(factor(condNames))[x]))
  condList <- lapply(indices, function(x) condNum[x])

  simpRES <- lapply(indices, function(x)
    alrInv(makeMat(x, model, bio = TRUE, useCov)))

  proportions <- t(mapply(function(x, y) x[[1]]$Estimate[match(uCond, c(refC, y))], simpRES, condList))
  colnames(proportions) <- paste0("Est_BioID", uCond)
  pVars <- t(mapply(function(x, y) x[[1]]$Variance[match(uCond, c(refC, y))], simpRES, condList))
  colnames(pVars) <- paste0("Var_BioID", uCond)
  pLL95 <- t(mapply(function(x, y) x[[1]]$LL95[match(uCond, c(refC, y))], simpRES, condList))
  colnames(pLL95) <- paste0("LL95_BioID", uCond)
  pUL95 <- t(mapply(function(x, y) x[[1]]$UL95[match(uCond, c(refC, y))], simpRES, condList))
  colnames(pUL95) <- paste0("UL95_BioID", uCond)

  bioDf <- data.frame(Protein = levels(factor(condNames)), proportions, pVars, pLL95, pUL95,
                       stringsAsFactors = F)

  #make avgCond log-ratio tables

  uBio <- levels(factor(oneDat$condID))
  condNum <- getCond(uBio, bio = FALSE)
  condNames <- getName(uBio)

  refC <- dat[[1]][1, "tag1"]

  nCond <- length(unique(condNum))
  uCond <- unique(c(refC, unique(condNum)))
  uCond <- uCond[order(uCond)]

  nProts <- length(unique(condNames))

  indices <- lapply(1:nProts, function(x)
    which(condNames == levels(factor(condNames))[x]))

  condList <- lapply(indices, function(x) condNum[x])

  simpRES <- lapply(indices, function(x)
    alrInv(makeMat(x, model, bio = FALSE, useCov, avgCond = TRUE)))


  refPos <- which(uCond == refC)
  suCond <- uCond[-refPos]
  if(length(suCond) == 1){
    lrs <- as.matrix(mapply(function(x, y) x[[2]]$Estimate[match(suCond, y)], simpRES, condList))
    colnames(lrs) <- paste0("Est_avg_fc", suCond)
    lrVars <- as.matrix(mapply(function(x, y) x[[2]]$Variance[match(suCond, y)], simpRES, condList))
    colnames(lrVars) <- paste0("Var_Cond", suCond)
    lrLL95 <- as.matrix(mapply(function(x, y) x[[2]]$LL95[match(suCond, y)], simpRES, condList))
    colnames(lrLL95) <- paste0("LL95_Cond", suCond)
    lrUL95 <- as.matrix(mapply(function(x, y) x[[2]]$UL95[match(suCond, y)], simpRES, condList))
    colnames(lrUL95) <- paste0("UL95_Cond", suCond)
  }else{
    lrs <- t(mapply(function(x, y) x[[2]]$Estimate[match(suCond, y)], simpRES, condList))
    colnames(lrs) <- paste0("Est_avg_fc", suCond)
    lrVars <- t(mapply(function(x, y) x[[2]]$Variance[match(suCond, y)], simpRES, condList))
    colnames(lrVars) <- paste0("Var_Cond", suCond)
    lrLL95 <- t(mapply(function(x, y) x[[2]]$LL95[match(suCond, y)], simpRES, condList))
    colnames(lrLL95) <- paste0("LL95_Cond", suCond)
    lrUL95 <- t(mapply(function(x, y) x[[2]]$UL95[match(suCond, y)], simpRES, condList))
    colnames(lrUL95) <- paste0("UL95_Cond", suCond)
  }
  avgLrTab <- data.frame(Protein = levels(factor(condNames)), lrs, lrVars, lrLL95, lrUL95,
                         stringsAsFactors = F)

  #Now get the proportions
  lrs <- t(mapply(function(x, y) x[[1]]$Estimate[match(uCond, unique(c(refC, y)))], simpRES, condList))
  colnames(lrs) <- paste0("Est_avg_Prop", uCond)
  lrVars <- t(mapply(function(x, y) x[[1]]$Variance[match(uCond, unique(c(refC, y)))], simpRES, condList))
  colnames(lrVars) <- paste0("Var_Cond", uCond)
  lrLL95 <- t(mapply(function(x, y) x[[1]]$LL95[match(uCond, unique(c(refC, y)))], simpRES, condList))
  colnames(lrLL95) <- paste0("LL95_Cond", uCond)
  lrUL95 <- t(mapply(function(x, y) x[[1]]$UL95[match(uCond, unique(c(refC, y)))], simpRES, condList))
  colnames(lrUL95) <- paste0("UL95_Cond", uCond)

  avgPTab <- data.frame(Protein = levels(factor(condNames)), lrs, lrVars, lrLL95, lrUL95,
                        stringsAsFactors = F)


    #redo with population level tables
  uBio <- levels(factor(oneDat$condID))
  condNum <- getCond(uBio, bio = FALSE)
  condNames <- getName(uBio)

  refC <- dat[[1]][1, "tag1"]

  nCond <- length(unique(condNum))
  uCond <- unique(c(refC, unique(condNum)))
  uCond <- uCond[order(uCond)]

  nProts <- length(unique(condNames))

  indices <- lapply(1:nProts, function(x)
    which(condNames == levels(factor(condNames))[x]))

  condList <- lapply(indices, function(x) condNum[x])

  simpRES <- lapply(indices, function(x)
    alrInv(makeMat(x, model, bio = FALSE, useCov, avgCond = FALSE)))


  refPos <- which(uCond == refC)
  suCond <- uCond[-refPos]
  if(length(suCond) == 1){
    lrs <- as.matrix(mapply(function(x, y) x[[2]]$Estimate[match(suCond, y)], simpRES, condList))
    colnames(lrs) <- paste0("Est_avg_fc", suCond)
    lrVars <- as.matrix(mapply(function(x, y) x[[2]]$Variance[match(suCond, y)], simpRES, condList))
    colnames(lrVars) <- paste0("Var_Cond", suCond)
    lrLL95 <- as.matrix(mapply(function(x, y) x[[2]]$LL95[match(suCond, y)], simpRES, condList))
    colnames(lrLL95) <- paste0("LL95_Cond", suCond)
    lrUL95 <- as.matrix(mapply(function(x, y) x[[2]]$UL95[match(suCond, y)], simpRES, condList))
    colnames(lrUL95) <- paste0("UL95_Cond", suCond)
  }else{
    lrs <- t(mapply(function(x, y) x[[2]]$Estimate[match(suCond, y)], simpRES, condList))
    colnames(lrs) <- paste0("Est_avg_fc", suCond)
    lrVars <- t(mapply(function(x, y) x[[2]]$Variance[match(suCond, y)], simpRES, condList))
    colnames(lrVars) <- paste0("Var_Cond", suCond)
    lrLL95 <- t(mapply(function(x, y) x[[2]]$LL95[match(suCond, y)], simpRES, condList))
    colnames(lrLL95) <- paste0("LL95_Cond", suCond)
    lrUL95 <- t(mapply(function(x, y) x[[2]]$UL95[match(suCond, y)], simpRES, condList))
    colnames(lrUL95) <- paste0("UL95_Cond", suCond)
  }
  popLrTab <- data.frame(Protein = levels(factor(condNames)), lrs, lrVars, lrLL95, lrUL95,
                         stringsAsFactors = F)

  #Now get the proportions
  lrs <- t(mapply(function(x, y) x[[1]]$Estimate[match(uCond, unique(c(refC, y)))], simpRES, condList))
  colnames(lrs) <- paste0("Est_avg_Prop", uCond)
  lrVars <- t(mapply(function(x, y) x[[1]]$Variance[match(uCond, unique(c(refC, y)))], simpRES, condList))
  colnames(lrVars) <- paste0("Var_Cond", uCond)
  lrLL95 <- t(mapply(function(x, y) x[[1]]$LL95[match(uCond, unique(c(refC, y)))], simpRES, condList))
  colnames(lrLL95) <- paste0("LL95_Cond", uCond)
  lrUL95 <- t(mapply(function(x, y) x[[1]]$UL95[match(uCond, unique(c(refC, y)))], simpRES, condList))
  colnames(lrUL95) <- paste0("UL95_Cond", uCond)

  popPTab <- data.frame(Protein = levels(factor(condNames)), lrs, lrVars, lrLL95, lrUL95,
                        stringsAsFactors = F)


  } #end else (not simple case)

  RES <- list()
  RES[[1]] <- model
  RES[[2]] <- bioDf
  RES[[3]] <- avgLrTab
  RES[[4]] <- avgPTab
  RES[[5]] <- popLrTab
  RES[[6]] <- popPTab

  #add Gene back in
  if(!is.null(avgLrTab)){
    tempTab <- avgLrTab
    }else{
      tempTab <- bioDf
    }
  Gene <- oneDat$gene[match(tempTab$Protein, oneDat$protein)]
  gRES <- lapply(RES, function(x) if(is.data.frame(x)){
    data.frame(Gene, x)}else{x})

  gRES
} #end of compBayes function




#' Changing the comparisons of interest
#'
#' This function manipulates a Stanfit object to change the contrasts that
#' are estimated.
#'
#' @export
#' @param res The results object to be altered.
#' @param contrastMat A matrix specifying the contrasts of interest.
#'   each row represents a contrast and the number of columns must match
#'   the number of parameters.  If this is not specified then by default
#'   tables for all possible pairwise comparisons will be generated.
#'
contrastEst <- function(res, contrastMat = NULL, useCov = FALSE){
  #How many conditions are there?
  lrIndex <- grep("Est", colnames(res[[3]]))
  condIndex <- grep("Est", colnames(res[[4]]))
  nCond <- length(condIndex)
  protNames <- res[[3]]$Protein
  nProts <- length(protNames)

  #Check to see if contrastMat has the correct dimensions
  if(!is.null(contrastMat)){
    if(ncol(contrastMat) != nCond){
      stop("Specified contrast dimensions are incorrect")
    }
  }

  #determine model structure

  #Is it a population level model?
  if(is.null(res[[2]])){
    simpleMod = TRUE
  }else{
      simpleMod = FALSE
    }



  lrCond <- as.integer(substring(colnames(res[[3]])[lrIndex],
                                     regexpr("fc" , colnames(res[[3]])[lrIndex]) + 2))
  fullCond <- as.integer(substring(colnames(res[[4]])[condIndex],
                                 regexpr("op" , colnames(res[[4]])[condIndex]) + 2))
  refCond <- setdiff(fullCond, lrCond)

  #Need to find the parameter index for each change
  startPos <- c(0, cumsum(apply(res[[3]][ , lrIndex], 1, function(x) sum(!is.na(x))))) + 1

  indices <- list()
  for(i in 1:(nProts)){
      indices[[i]] <- startPos[i]:(startPos[i + 1] - 1)
  }

  #Find set of observed conditions for each protein
  condList <- lapply(1:nProts, function(x)
    fullCond[which(!is.na(res[[4]])[x, condIndex])])

  #make simplexes
  if(simpleMod){
    avgSimp <- lapply(indices, function(x)
      alrInv(makeMat(x, res[[1]], bio = FALSE, useCov, avgCond = FALSE), justSimp = TRUE))
  }else{
    avgSimp <- lapply(indices, function(x)
      alrInv(makeMat(x, res[[1]], bio = FALSE, useCov, avgCond = TRUE), justSimp = TRUE))

    popSimp <- lapply(indices, function(x)
      alrInv(makeMat(x, res[[1]], bio = FALSE, useCov, avgCond = FALSE), justSimp = TRUE))
  }

  newRef <- list()
  for(j in 1:(length(fullCond) - 1)){
    newList <- lapply(1:nProts, function(x)
      summSimp(avgSimp[[x]], fullCond, condList[[x]], lrCond[j], protNames[x]))
    avgLrTab <- do.call(rbind, newList)

    if(simpleMod == FALSE){
      newList <- lapply(1:nProts, function(x)
        summSimp(popSimp[[x]], fullCond, condList[[x]], lrCond[j], protNames[x]))
      popLrTab <- do.call(rbind, newList)
    }else{
      popLrTab <- NULL
    }

    newRef[[j]] <- list(avgLrTab, popLrTab)
  }

#all pairwise comparisons are done.  Now make contrasts if any were specified.

if(is.null(contrastMat)){
  contRes <- NULL
}else{
  nCont <- nrow(contrastMat)
  contRes <- list()
  for(k in 1:nCont){
    tempList <- lapply(1:nProts, function(x)
      contSimp(avgSimp[[x]], fullCond, condList[[x]], contrastMat[k, ], protNames[x]))
    avgCont <- do.call(rbind, tempList)


    if(simpleMod == FALSE){
      tempList <- lapply(1:nProts, function(x)
        contSimp(popSimp[[x]], fullCond, condList[[x]], contrastMat[k, ], protNames[x]))
      popCont <- do.call(rbind, tempList)
    }else{
      popCont <- NULL
    }

    contRes[[k]] <- list(avgCont, popCont)
  }

}

list(newRef, contRes)

}#end contrastEst function



