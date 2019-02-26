#Functions to help with longitudinal analysis

#' This function runs the main code for fitting a model that compares more
#' than one curve through time.
#'
#' @export
#' @param tempDat The data to be analyzed
#' @param timeDegree either 1, 2 or 3 for linear, quadratic or cubic models
#' @param fullTimes A vector specifying all of the times found in the study
#' @param useW A boolean that determines whether or not the posterior variance of
#'   each protein estimate should be used to weight the observations
#' @param refCat A string that specifies the reference category
#' @param groupByGene A boolean that specifies whether timepoints are grouped
#'   by genes or proteins
#'
testInteract <- function(tempDat, timeDegree = 2, fullTimes, useW = TRUE,
                  refCat = NULL, groupByGene = FALSE){


  if(groupByGene){
    uProt <- unique(tempDat$Gene)
    #Model always uses Protein column
    #Replace column and save the protein names for later
    savedProts <- tempDat$Protein
    tempDat$Protein <- tempDat$Gene
  }else{
    uProt <- unique(tempDat$Protein)
    savedProts <- tempDat$Protein
  }

  nProt <- length(uProt)
  obsTimes <- unique(tempDat$Time)


  degVec <- as.character(1:timeDegree)
  degVec[1] <- ""

  #fmla <- as.formula(paste("FC ~ 0 + Gene + Gene:Species + Gene:Time + Gene:Time2 +
   #                        Gene:Species:Time + Gene:Species:Time2 "))

  fmla <- as.formula(paste("FC ~ 0 + Protein + Protein:Category + ",
                           paste0("Protein:", paste0("Time", degVec), collapse = " + "),
                            "+", paste0("Protein:Category:", paste0("Time", degVec), collapse = " + ")))

  if(!is.null(refCat)){
    tempDat$Category <- factor(tempDat$Category)
    tempDat$Category <- relevel(tempDat$Category, ref = refCat)
  }else{
    tempDat$Category <- factor(tempDat$Category)
  }

  uCats <- levels(tempDat$Category)

  #Create strings that define the model
  if(useW){
    fullMod <- lm(fmla, weights = w_, data = tempDat)
  }else{
    fullMod <- lm(fmla, data = tempDat)
  }


  #create list of times within observed ranges
  times <- list()
  for(c_ in 1:length(uCats)){
    catTime <- unique(tempDat[which(tempDat$Category == uCats[c_]), "Time"])
    maxT <- max(catTime)
    minT <- min(catTime)
    boolT <- (obsTimes >= minT) & (obsTimes <= maxT)
    smalltime <- obsTimes[boolT]
    times[[c_]] <- smalltime
  }

  #Now get predictions and p-values for each protein
  pRes <- list()

  for(index in 1:nProt){
    #Implement F test for time effect
    #Create strings that define the hypothesis tests



    timeTests <- list()
    catTests <- list()

    #Test for an overall time effect in the baseline condition
    timeStr <-  paste0("Protein", uProt[index], paste0(":Time", degVec), " = 0")
    timeTests[[1]] <-  lht(fullMod, timeStr)$`Pr(>F)`[2]


    for(t_ in 2:length(uCats)){
      #Test for any difference between refCat and category t_
      catStr <-  paste0("Protein", uProt[index], ":Category", uCats[t_], c("", paste0(":Time", degVec)), " = 0")
      catTests[[t_ - 1]] <- lht(fullMod, catStr)$`Pr(>F)`[2]

      #Test for an overall time effect in condition t_
      timeStr <- paste0("Protein", uProt[index], paste0(":Time", degVec), " + ",
                        "Protein", uProt[index], ":Category", uCats[t_], paste0(":Time", degVec)," = 0")
      timeTests[[t_]] <- lht(fullMod, timeStr)$`Pr(>F)`[2]
    }

    #Now extract and store the predicted values
    newDfs <- list()
    if(timeDegree == 1){
      newDfs <- lapply(1:length(uCats), function(x) data.frame(Protein = uProt[index], Time = times[[x]], Category = uCats[x]))
    }
    if(timeDegree == 2){
      newDfs <- lapply(1:length(uCats), function(x) data.frame(Protein = uProt[index], Time = times[[x]], Time2 = times[[x]]^2,
                                                               Category = uCats[x]))
    }
    if(timeDegree == 3){
      newDfs <- lapply(1:length(uCats), function(x) data.frame(Protein = uProt[index], Time = times[[x]], Time2 = times[[x]]^2,
                                                               Time3 = times[[x]]^3, Category = uCats[x]))
    }

    catPreds <- lapply(newDfs, function(x) predict(fullMod, x))

    catRes <- list()
    catRes[[1]] <- data.frame(Gene = tempDat$Gene[match(uProt[index], tempDat$Protein)],
                              Protein = savedProts[match(uProt[index], tempDat$Protein)],
                              matrix(catPreds[[1]][match(fullTimes, times[[1]])], nrow = 1), timeTests[[1]])
    colnames(catRes[[1]]) <- c("Gene", "Protein", paste("category:",uCats[1], paste0("Time", fullTimes)), paste("category:",uCats[1], "Pval-Time"))

    for(t_ in 2:length(uCats)){
      catRes[[t_]] <- data.frame(matrix(catPreds[[t_]][match(fullTimes, times[[t_]])], nrow = 1), timeTests[[t_]], catTests[[t_ - 1]])
      colnames(catRes[[t_]]) <- c(paste("category:",uCats[t_], paste0("Time", fullTimes)), paste("category:", uCats[t_], "Pval-Time"), paste0("Pval-", uCats[t_]))
    }

    pRes[[index]] <- do.call(cbind, catRes)

  }#end Gene loop

  resDf <- do.call(rbind, pRes)

  #Now add columns with q values
  #first find the p-values

  pIndex <- grep("Pval", colnames(resDf))
  Qvals <- as.data.frame(lapply(resDf[ , pIndex], function(x) p.adjust(x, method = "fdr")))
  colnames(Qvals) <- gsub("Pval", "Qval", colnames(resDf)[pIndex])
  finalDf <- data.frame(resDf, Qvals)

  newIndex <- c(seq_along(resDf), pIndex + 0.5)
  finalDf <- finalDf[ ,order(newIndex)]

  finalDf
}

#' This function runs the main code for fitting a model that looks at the
#' overall (single curve) effect through time.
#'
#' @export
#' @param tempDat The data to be analyzed
#' @param timeDegree either 1, 2 or 3 for linear, quadratic or cubic models
#' @param fullTimes A vector specifying all of the times found in the study
#' @param useW A boolean that determines whether or not the posterior variance of
#'   each protein estimate should be used to weight the observations
#'
# OVERALL EFFECT
test_overall_effect <- function(tempDat, timeDegree = 2, fullTimes, useW = TRUE){

  uProt <- unique(tempDat$Protein)
  nProt <- length(uProt)
  times <- unique(tempDat$Time)

  degVec <- as.character(1:timeDegree)
  degVec[1] <- ""

  fmla <- as.formula(paste("FC ~ Protein * (",
                     paste0("Time", degVec, collapse = " + "), ")",
                         " - ", paste0("Time", degVec, collapse = " - ")))



  if(useW){
    fullMod <- lm(fmla, weights = w_, data = tempDat)
  }else{
    fullMod <- lm(fmla, data = tempDat)
  }

  # Now do a test for time effects on each protein
  pList <- list()

  for(index in 1:nProt){
    # Implement F test for time effect
    testl <- lht(fullMod, paste(paste0("Protein", uProt[index], ":", paste0("Time", degVec), " = 0")))
    pVal <- testl$'Pr(>F)'[2]
    # Now extract and store the predicted values
    if(timeDegree == 1){
      newDf <- data.frame(Protein = uProt[index], Time = times)
    }
    if(timeDegree == 2){
      newDf <- data.frame(Protein = uProt[index], Time = times, Time2 = times^2)
    }
    if(timeDegree == 3){
      newDf <- data.frame(Protein = uProt[index], Time = times, Time2 = times^2, Time3 = times^3)
    }

    mFit <- predict(fullMod, newDf)
    protList <- list(as.character(uProt[index]), c(mFit[match(fullTimes, times)], pVal))
    names(protList[[2]]) <- c(paste0("Time", fullTimes), "Pval")
    pList[[index]] <- protList
  } # end protein loop

  resNames <- unlist(lapply(pList, function(x) x[[1]]))
  resVals <- do.call(rbind, lapply(pList, function(x) x[[2]]))
  resDf <- data.frame(Protein = resNames, resVals)

  pIndex <- grep("Pval", colnames(resDf))
  Qvals <- p.adjust(resDf[ , pIndex], method = "fdr")

  resDf$Qval <- Qvals
  gRes <- data.frame(Gene = tempDat$Gene[match(resDf$Protein, tempDat$Protein)], resDf)

  gRes
}
