#File for manipulating the results of a model with PTMs

#function for extracting alpha LRs and Props

getPtms <- function(model, ptmDat, refC){

  uBio <- levels(factor(ptmDat$ptmID))
  condNum <- getCond(uBio, bio = FALSE)
  protNames <- getName(uBio)
  Gene <- ptmDat$gene[match(protNames, ptmDat$protein)]

  # ensures peptide names stay the same regardless of condition number size
  Peptide <- unlist(lapply(uBio, function(x) paste(strsplit(x, "_")[[1]][3], strsplit(x, "_")[[1]][4], sep="_")))

  nCond <- length(unique(condNum))
  uCond <- unique(c(refC, unique(condNum)))
  uCond <- uCond[order(uCond)]

  nPeps <- length(unique(Peptide))
  uPep <- levels(factor(Peptide))

  indices <- lapply(1:nPeps, function(x)
    which(Peptide == levels(factor(Peptide))[x]))

  condList <- lapply(indices, function(x) condNum[x])

  simpRES <- lapply(indices, function(x) {
    extracted <- lapply(x, function(y)
      rstan::extract(model,
                     pars=paste("alpha[", y, "]", sep = ""))$alpha)
    extractMat <- do.call(cbind, extracted)
    alrInv(extractMat)
    })


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
  avgLrTab <- data.frame(Gene = Gene[match(uPep, Peptide)], Protein = protNames[match(uPep, Peptide)],
                         Peptide = uPep, lrs, lrVars, lrLL95, lrUL95,
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

  avgPTab <- data.frame(Gene = Gene[match(uPep, Peptide)], Protein = protNames[match(uPep, Peptide)],
                        Peptide = uPep, lrs, lrVars, lrLL95, lrUL95,
                        stringsAsFactors = F)

list(avgLrTab, avgPTab)

}




#A function for extracting all contrasts from the PTM parameters
contrastPtm <- function(res, contrastMat = NULL){
  #How many conditions are there?
  lrIndex <- grep("Est", colnames(res[[9]]))
  condIndex <- grep("Est", colnames(res[[10]]))
  nCond <- length(condIndex)
  # If there's only 2 conditions, no extra comparisons can be made, no need to proceed
  if(nCond == 2){
    stop("Only 2 conditions detected, no extra comparisons to be made.")
  }

  pepNames <- res[[9]]$Peptide

  nPeps <- length(pepNames)

  #Check to see if contrastMat has the correct dimensions
  if(!is.null(contrastMat)){
    if(ncol(contrastMat) != nCond){
      stop("Specified contrast dimensions are incorrect")
    }
  }

  #read in relevant variables
  ptmDat <- res[[8]]
  model <- res[[1]]

  #follow the same code used for the original simplex creation (but use simpOnly)
  lrCond <- as.integer(substring(colnames(res[[9]])[lrIndex],
                                 regexpr("Fc" , colnames(res[[9]])[lrIndex]) + 2))
  fullCond <- as.integer(substring(colnames(res[[10]])[condIndex],
                                   regexpr("op" , colnames(res[[10]])[condIndex]) + 2))
  refCond <- setdiff(fullCond, lrCond)

  uBio <- levels(factor(ptmDat$ptmID))
  condNum <- getCond(uBio, bio = FALSE)
  protNames <- getName(uBio)
  Gene <- ptmDat$gene[match(protNames, ptmDat$protein)]

  Peptide <- unlist(lapply(uBio, function(x) paste(strsplit(x, "_")[[1]][3], strsplit(x, "_")[[1]][4], sep="_")))

  nCond <- length(unique(condNum))
  uCond <- unique(c(refCond, unique(condNum)))
  uCond <- uCond[order(uCond)]

  nPeps <- length(unique(Peptide))
  uPep <- levels(factor(Peptide))

  indices <- lapply(1:nPeps, function(x)
    which(Peptide == levels(factor(Peptide))[x]))

  condList <- lapply(indices, function(x) condNum[x])

  simpRES <- lapply(indices, function(x) {
    extracted <- lapply(x, function(y)
      rstan::extract(model,
                     pars=paste("alpha[", y, "]", sep = ""))$alpha)
    extractMat <- do.call(cbind, extracted)
    alrInv(extractMat, justSimp = TRUE)
  })


  newRef <- list()
  for(j in 1:(length(fullCond) - 1)){
    newList <- lapply(1:nPeps, function(x)
      summSimp(simpRES[[x]], fullCond, condList[[x]], lrCond[j], uPep[x]))
    avgLrTab <- do.call(rbind, newList)
    colnames(avgLrTab)[1] <- "Peptide"
    avgLrTab <- data.frame(Gene = res[[9]]$Gene, Protein = res[[9]]$Protein, avgLrTab)

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
      tempList <- lapply(1:nPeps, function(x)
        contSimp(avgSimp[[x]], fullCond, condList[[x]], contrastMat[k, ], uPep[x]))
      avgCont <- do.call(rbind, tempList)
      colnames(avgCont)[1] <- "Peptide"
      avgCont <- data.frame(Gene = res[[9]]$Gene, Protein = res[[9]]$Protein, avgCont)


      popCont <- NULL


      contRes[[k]] <- list(avgCont, popCont)
    }

  }

  list(newRef, contRes)

}#end contrastPtm function
