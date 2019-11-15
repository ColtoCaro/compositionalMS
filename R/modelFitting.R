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
#' @param bridge If bridge == TRUE then the first column will be treated as the
#'   reference for the entire experiment.  If the condition number repeats, these
#'   entries will be collapsed into a single reference channel.  When bridge is false
#'   this collapse will prevent estimation of the individual collapsed channels which
#'   are necessary for viewing proportion plots.  In this case a second model will be
#'   fit for the sole purpose of visualizing the behavior of individual replicates.
#' @param adapt_delta The target proposal acceptance rate that determines step
#'   size in a Stan model.
#' @param varPool A variable that determines the structure of the variance parameters.
#'   varPool = 0 denotes no pooling, i.e. a separate variance parameter for each protein.
#'   varPool = 1 (default) generates partially pooled variance parameters and
#'   varPool = 2 creates complete pooling, i.e. only one experimental error.
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
                     normalize_ptm = FALSE,
                     pop_sd = 10,
                     simpleMod = FALSE,
                     bridge = TRUE,
                     adapt_delta = .9,
                     varPool = 1,
                     no_bio = FALSE                     ){

  #Force the new setup with two simple runs
  simpleMod <- TRUE
  bridge <- FALSE

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


  #setup two runs of the function if conds == bioreps
  condVec <- unlist(lapply(dat, function(x) x[1, grep("tag", colnames(x))]))
  bioVec <- unlist(lapply(dat, function(x) x[2, grep("tag", colnames(x))]))

  if(sum(condVec == bioVec) < length(condVec)){
    model_number <- 2
    simpleMod <- TRUE
  }else{
      model_number <- 1
    }

  for(modelN in 1:model_number){

  #If this is the second time through, set bioReps to conditions
  if(modelN == 2){
    dat <- lapply(dat, function(x){
      x[1, ] <- x[2, ]
      x
      })
    ptmTemp <- oneDat$ptmID #save this from the first model
  }
  #how to normalize
  if (normalize) {
    normalize_plex <- lapply(1:length(dat), function(x) (all(dat[[x]][3,] == 0)))
    }
  if (normalize & normalize_ptm) {
    normalize_plex <- lapply(1:length(dat), function(x) (TRUE))
    }
  if ( ! (normalize | normalize_ptm)){
    normalize_plex <- lapply(1:length(dat), function(x) (FALSE))
    }
  readyDat <- lapply(1:length(dat), function(x)
    transformDat(dat[[x]], plexNumber = x, normalize = normalize_plex[[x]],
                 simpleMod))
  oneDat <- do.call(rbind, readyDat)
  if ( !  normalize_ptm){
  oneDat <- oneDat[oneDat$lr != 0,]
  }
  #set data variables
  #Do ptms first since it might change the dimensions of the data
  #N_ <- nrow(oneDat)
  #count the plexes with PTMs
  sumPtm <- sum(unlist(lapply(dat, function(x) (sum(x[3, grep("tag", colnames(x))]) > 0))))

  #If this is the second time through, set bioReps to conditions
  #First time through, remove PTMs
  if(modelN ==2){
    oneDat$condID <- oneDat$bioID
  } #else{
    #if(sumPtm > 0){oneDat <- oneDat[which(oneDat$ptm == 0), ] }
  #}

  N_ <- nrow(oneDat)
  if(sumPtm == 0){
    n_p <- 0
    n_ptm <- 0
    ptm <- rep(0,N_)
    ptmPep <- rep(0,N_)
  }else{
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
      #temporary reduction for estimating only PTMs
        globOnly <- setdiff(globalProts, ptmProts)
        globIndex <- which(oneDat$bioID %in% globOnly)
        if (length(globIndex) > 0){oneDat <- oneDat[-globIndex, ]}
    }

    nonPtms <- which(oneDat$ptm == 0)
    ptmDat <- oneDat[-nonPtms, ]
    ptmName <- levels(factor(ptmDat$ptmID))
    n_p <- length(ptmName)

    n_ptm <- length(unique(ptmDat$ptm))
    # if you want variance estimate per peptide per PTM change to the line below
    # n_ptm <- length(unique(ptmDat$ptmID))

    ptm <- as.integer(oneDat$ptm)
    ptmPep <- rep(0, nrow(oneDat))
    ptmPep[which(oneDat$ptm > 0)] <- as.integer(factor(ptmDat$ptmID))

  } # end actions for ptm experiments

  oneDat <- oneDat[order(oneDat$condID, oneDat$bioID, oneDat$ptm,
                         oneDat$ptmID), ]
  # remove log
  # start setting global variable for stan model
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
    n_pep <- table(condID)
  }else{
    n_b <- length(unique(oneDat$bioID))
    bioID <- as.integer(factor(oneDat$bioID))
    n_pep <- table(bioID)
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
  #model <- stan(file="exec/allModels.stan",
   #             iter = 2000, cores = 4, control = list(adapt_delta = .9))

  sMod <- compMS:::stanmodels$allModels
  model <- rstan::sampling(sMod, cores = nCores, iter = iter,
                           control = list(adapt_delta = adapt_delta))


  #create summary tables

  #simple case
  # if number of biological replicates equals number of conditions
  # or if simple model was set to TRUE
  if((n_b == n_c) | simpleMod){

      uBio <- levels(factor(oneDat$condID))
      if(modelN == 1){
        condNum <- getCond(uBio, bio = FALSE)
      }else{
        condNum <- getCond(uBio, bio = TRUE)
      }
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
        colnames(lrs) <- paste0("Est_Fc", suCond)
        lrVars <- as.matrix(mapply(function(x, y) x[[2]]$Variance[match(suCond, y)], simpRES, condList))
        colnames(lrVars) <- paste0("Var_", suCond)
        lrLL95 <- as.matrix(mapply(function(x, y) x[[2]]$LL95[match(suCond, y)], simpRES, condList))
        colnames(lrLL95) <- paste0("LL95_", suCond)
        lrUL95 <- as.matrix(mapply(function(x, y) x[[2]]$UL95[match(suCond, y)], simpRES, condList))
        colnames(lrUL95) <- paste0("UL95_", suCond)
      }else{
        lrs <- t(mapply(function(x, y) x[[2]]$Estimate[match(suCond, y)], simpRES, condList))
        colnames(lrs) <- paste0("Est_Fc", suCond)
        lrVars <- t(mapply(function(x, y) x[[2]]$Variance[match(suCond, y)], simpRES, condList))
        colnames(lrVars) <- paste0("Var_", suCond)
        lrLL95 <- t(mapply(function(x, y) x[[2]]$LL95[match(suCond, y)], simpRES, condList))
        colnames(lrLL95) <- paste0("LL95_", suCond)
        lrUL95 <- t(mapply(function(x, y) x[[2]]$UL95[match(suCond, y)], simpRES, condList))
        colnames(lrUL95) <- paste0("UL95_", suCond)
      }
      avgLrTab <- data.frame(Protein = levels(factor(condNames)), lrs, lrVars, lrLL95, lrUL95,
                             stringsAsFactors = F)

      #Now get the proportions
      lrs <- t(mapply(function(x, y) x[[1]]$Estimate[match(uCond, unique(c(refC, y)))], simpRES, condList))
      colnames(lrs) <- paste0("Est_Prop", uCond)
      lrVars <- t(mapply(function(x, y) x[[1]]$Variance[match(uCond, unique(c(refC, y)))], simpRES, condList))
      colnames(lrVars) <- paste0("Var_", uCond)
      lrLL95 <- t(mapply(function(x, y) x[[1]]$LL95[match(uCond, unique(c(refC, y)))], simpRES, condList))
      colnames(lrLL95) <- paste0("LL95_", uCond)
      lrUL95 <- t(mapply(function(x, y) x[[1]]$UL95[match(uCond, unique(c(refC, y)))], simpRES, condList))
      colnames(lrUL95) <- paste0("UL95_", uCond)

      avgPTab <- data.frame(Protein = levels(factor(condNames)), lrs, lrVars, lrLL95, lrUL95,
                            stringsAsFactors = F)


      bioDf <- NULL
      fcBioTab <- NULL
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
  colnames(proportions) <- paste0("Est_Prop", uCond)
  pVars <- t(mapply(function(x, y) x[[1]]$Variance[match(uCond, c(refC, y))], simpRES, condList))
  colnames(pVars) <- paste0("Var_", uCond)
  pLL95 <- t(mapply(function(x, y) x[[1]]$LL95[match(uCond, c(refC, y))], simpRES, condList))
  colnames(pLL95) <- paste0("LL95_", uCond)
  pUL95 <- t(mapply(function(x, y) x[[1]]$UL95[match(uCond, c(refC, y))], simpRES, condList))
  colnames(pUL95) <- paste0("UL95_", uCond)

  bioDf <- data.frame(Protein = levels(factor(condNames)), proportions, pVars, pLL95, pUL95,
                       stringsAsFactors = F)

  #Now do fcBioTab
  uCond <- unique(condNum)
  uCond <- uCond[order(uCond)]

  lrs <- t(mapply(function(x, y) x[[2]]$Estimate[match(uCond, y)], simpRES, condList))
  colnames(lrs) <- paste0("Est_Fc", uCond)
  pVars <- t(mapply(function(x, y) x[[2]]$Variance[match(uCond, y)], simpRES, condList))
  colnames(pVars) <- paste0("Var_", uCond)
  pLL95 <- t(mapply(function(x, y) x[[2]]$LL95[match(uCond, y)], simpRES, condList))
  colnames(pLL95) <- paste0("LL95_", uCond)
  pUL95 <- t(mapply(function(x, y) x[[2]]$UL95[match(uCond, y)], simpRES, condList))
  colnames(pUL95) <- paste0("UL95_", uCond)

  fcBioTab <- data.frame(Protein = levels(factor(condNames)), lrs, pVars, pLL95, pUL95,
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
    colnames(lrs) <- paste0("Est_Fc", suCond)
    lrVars <- as.matrix(mapply(function(x, y) x[[2]]$Variance[match(suCond, y)], simpRES, condList))
    colnames(lrVars) <- paste0("Var_", suCond)
    lrLL95 <- as.matrix(mapply(function(x, y) x[[2]]$LL95[match(suCond, y)], simpRES, condList))
    colnames(lrLL95) <- paste0("LL95_", suCond)
    lrUL95 <- as.matrix(mapply(function(x, y) x[[2]]$UL95[match(suCond, y)], simpRES, condList))
    colnames(lrUL95) <- paste0("UL95_", suCond)
  }else{
    lrs <- t(mapply(function(x, y) x[[2]]$Estimate[match(suCond, y)], simpRES, condList))
    colnames(lrs) <- paste0("Est_Fc", suCond)
    lrVars <- t(mapply(function(x, y) x[[2]]$Variance[match(suCond, y)], simpRES, condList))
    colnames(lrVars) <- paste0("Var_", suCond)
    lrLL95 <- t(mapply(function(x, y) x[[2]]$LL95[match(suCond, y)], simpRES, condList))
    colnames(lrLL95) <- paste0("LL95_", suCond)
    lrUL95 <- t(mapply(function(x, y) x[[2]]$UL95[match(suCond, y)], simpRES, condList))
    colnames(lrUL95) <- paste0("UL95_", suCond)
  }
  avgLrTab <- data.frame(Protein = levels(factor(condNames)), lrs, lrVars, lrLL95, lrUL95,
                         stringsAsFactors = F)

  #Now get the proportions
  lrs <- t(mapply(function(x, y) x[[1]]$Estimate[match(uCond, unique(c(refC, y)))], simpRES, condList))
  colnames(lrs) <- paste0("Est_Prop", uCond)
  lrVars <- t(mapply(function(x, y) x[[1]]$Variance[match(uCond, unique(c(refC, y)))], simpRES, condList))
  colnames(lrVars) <- paste0("Var_", uCond)
  lrLL95 <- t(mapply(function(x, y) x[[1]]$LL95[match(uCond, unique(c(refC, y)))], simpRES, condList))
  colnames(lrLL95) <- paste0("LL95_", uCond)
  lrUL95 <- t(mapply(function(x, y) x[[1]]$UL95[match(uCond, unique(c(refC, y)))], simpRES, condList))
  colnames(lrUL95) <- paste0("UL95_", uCond)

  avgPTab <- data.frame(Protein = levels(factor(condNames)), lrs, lrVars, lrLL95, lrUL95,
                        stringsAsFactors = F)


  } #end else (not simple case)

  if(sumPtm > 0){
    ptmTabs <- getPtms(model, ptmDat, dat[[1]][1, "tag1"])
    ptmLr <- ptmTabs[[1]]
    ptmProp <- ptmTabs[[2]]
  }else{
    ptmDat <- NULL
    ptmLr <- NULL
    ptmProp <- NULL
  }


  if(modelN == 1){
    RES <- list()

    RES[[1]] <- model
    RES[[2]] <- bioDf
    RES[[3]] <- avgLrTab
    RES[[4]] <- avgPTab
    RES[[5]] <- fcBioTab
    RES[[6]] <- NULL  #This used to be for population models
    RES[[7]] <- oneDat$condID
    RES[[8]] <- ptmDat
    RES[[9]] <- ptmLr
    RES[[10]] <- ptmProp
    RES[[11]] <- oneDat #useful for troubleshooting
  }else{
    RES[[2]] <- avgPTab
    RES[[5]] <- avgLrTab
    RES[[12]] <- model

  }

if (no_bio){
break
}
  }  #End modelN for-loop that fits 2 simple models

  #add Gene back in
  if(!is.null(avgLrTab)){
    tempTab <- avgLrTab
    }else{
      tempTab <- bioDf
    }
  # add genes column without worrying about gene columns being lost
  gRES <- RES
  x <- oneDat[match(tempTab$Protein, oneDat$protein), c('gene', 'protein')]
  for(i in 2:5){
    if(is.data.frame(RES[[i]])){
      y <- gRES[[i]]
      gRES[[i]] <- merge(x, y, by.x='protein', by.y='Protein', all=T)
      colnames(avgLrTab)[colnames(avgLrTab) == c('gene')] <-  'Gene'
      colnames(avgLrTab)[colnames(avgLrTab) == c('protein')] <- 'Protein'
    }
  }

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
contrastEst <- function(res, contrastMat = NULL, useCov = FALSE, stanKey = NULL){
  #How many conditions are there?
  lrIndex <- grep("Est", colnames(res[[3]]))
  condIndex <- grep("Est", colnames(res[[4]]))
  nCond <- length(condIndex)
  # If there's only 2 conditions, no extra comparisons can be made, no need to proceed
  if(nCond == 2){
    stop("Only 2 conditions detected, no extra comparisons to be made.")
  }

  protNames <- res[[3]]$Protein
  nProts <- length(protNames)

  #Check to see if contrastMat has the correct dimensions
  if(!is.null(contrastMat)){
    if(ncol(contrastMat) != nCond){
      stop("Specified contrast dimensions are incorrect")
    }
  }

  #determine model structure

  #Is it a population level model? -  Answer is now NO
  if(is.null(res[[2]])){
    simpleMod = TRUE
  }else{
      simpleMod = TRUE
    }


  lrCond <- as.integer(substring(colnames(res[[3]])[lrIndex],
                                     regexpr("Fc" , colnames(res[[3]])[lrIndex]) + 2))
  fullCond <- as.integer(substring(colnames(res[[4]])[condIndex],
                                 regexpr("op" , colnames(res[[4]])[condIndex]) + 2))
  refCond <- setdiff(fullCond, lrCond)

  #Need to find the parameter index for each change
  startPos <- c(0, cumsum(apply(res[[3]][ , lrIndex], 1, function(x) sum(!is.na(x))))) + 1

  if(is.null(stanKey)){
   stop("Please enter a key to match proteins with Stan variables")
  indices <- list()
    for(i in 1:(nProts)){
      indices[[i]] <- startPos[i]:(startPos[i + 1] - 1)
    }
  }else{
    uBio <- levels(factor(stanKey))
    condNum <- getCond(uBio, bio = FALSE)
    condNames <- getName(uBio)

    indices <- lapply(1:nProts, function(x)
      which(condNames == levels(factor(condNames))[x]))

    condList <- lapply(indices, function(x) condNum[x])
  }

  #Find set of observed conditions for each protein
  #condList <- lapply(1:nProts, function(x)
  #  fullCond[which(!is.na(res[[4]])[x, condIndex])])

  #make simplexes
  if(simpleMod){
    avgSimp <- lapply(indices, function(x)
      alrInv(makeMat(x, res[[1]], bio = FALSE, useCov, avgCond = FALSE), justSimp = TRUE))
  }else{
    avgSimp <- lapply(indices, function(x)
      alrInv(makeMat(x, res[[1]], bio = FALSE, useCov, avgCond = TRUE), justSimp = TRUE))

    #popSimp <- lapply(indices, function(x)
    #  alrInv(makeMat(x, res[[1]], bio = FALSE, useCov, avgCond = FALSE), justSimp = TRUE))
  }

  newRef <- list()
  for(j in 1:(length(fullCond) - 1)){
    newList <- lapply(1:nProts, function(x)
      summSimp(avgSimp[[x]], fullCond, condList[[x]], lrCond[j], protNames[x]))
    avgLrTab <- do.call(rbind, newList)
    avgLrTab <- data.frame(Gene = res[[3]]$Gene, avgLrTab)

    #if(simpleMod == FALSE){
     # newList <- lapply(1:nProts, function(x)
      #  summSimp(popSimp[[x]], fullCond, condList[[x]], lrCond[j], protNames[x]))
    #   popLrTab <- do.call(rbind, newList)
    #   popLrTab <- data.frame(Gene = res[[3]]$Gene, popLrTab)
    # }else{
      popLrTab <- NULL
    #}

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
    avgCont <- data.frame(Gene = res[[3]]$Gene, avgCont)

    # if(simpleMod == FALSE){
    #   tempList <- lapply(1:nProts, function(x)
    #     contSimp(popSimp[[x]], fullCond, condList[[x]], contrastMat[k, ], protNames[x]))
    #   popCont <- do.call(rbind, tempList)
    #   popCont <- data.frame(Gene = res[[3]]$Gene, popCont)
    # }else{
      popCont <- NULL
    #}

    contRes[[k]] <- list(avgCont, popCont)
  }

}

list(newRef, contRes)

}#end contrastEst function



