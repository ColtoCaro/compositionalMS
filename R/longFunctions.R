#Functions to help with longitudinal analysis

#' This function runs the main code for fitting a model that compares more
#' than one curve through time.
#'
#' @export
#' @param tempDat The data to be analyzed
#' @param timeDegree either 1, 2 or 3 for linear, quadratic or cubic models
#' @param fullTimes A vector specifying all of the times found in the study
#' @param fullCats A vector specifying all of the categories found in the study
#' @param useW A boolean that determines whether or not the posterior variance of
#'   each protein estimate should be used to weight the observations
#' @param refCat A string that specifies the reference category
#' @param groupByGene A boolean that specifies whether timepoints are grouped
#'   by genes or proteins
#' @param randEffect takes 0 for no random effects. 1 for a random intercept
#'   and 2 for a random slope and intercepts model
#'
testInteract <- function(tempDat, timeDegree = 2, fullTimes, fullCats, useW = TRUE,
                  refCat = NULL, groupByGene = FALSE, randEffect = 0){

  #Establish relevant grouping and make sure data is ordered
  if(groupByGene){
    #Model always uses Protein column
    tempDat <- tempDat[order(tempDat$Gene), ]
    #Save the protein names for later, and replace column
    savedProts <- tempDat$Protein
    tempDat$Protein <- tempDat$Gene
    uProt <- unique(tempDat$Gene)
  }else{

    tempDat <- tempDat[order(tempDat$Protein), ]
    savedProts <- tempDat$Protein
    uProt <- unique(tempDat$Protein)
  }

  #extract dimensions of the data
  nProt <- length(uProt)
  obsTimes <- unique(tempDat$Time)



  #Determine model details and adjust entries based on the column names
  randIndex <- grep("Random_Effect", colnames(tempDat))
  if(length(randIndex > 0)){
    if(length(randIndex) > 1){stop("Only one random intercept is allowed")}
    mixedMod <- TRUE
    #concatenate protein and random ID
    tempDat$Random <- paste0(tempDat$Protein, tempDat[ , randIndex])
  }else{
    mixedMod <- FALSE
  }

  contIndex <- grep("Continuous_Covariate", colnames(tempDat))
  if(length(contIndex > 0)){
    contCovar <- TRUE
    #center the continuous variables
    #for(i in 1:length(contIndex)){
    #  varMean <- mean(tempDat[ , contIndex[i]], na.rm = TRUE)
    #  tempDat[ , contIndex[i]] <- tempDat[ , contIndex[i]] - varMean
    #}
  }else{
    contCovar <- FALSE
  }

  catCovarIndex <- grep("Categorical_Covariate", colnames(tempDat))
  if(length(catCovarIndex > 0)){
    catCovar <- TRUE
    #make vector of references
    refList <- list()
    levelList <- list()
    #make sure these variables are treated as factors
    for(i in 1:length(catCovarIndex)){
      tempDat[ , catCovarIndex[i]] <- factor(tempDat[ , catCovarIndex[i]])
      refList[[i]] <- levels(tempDat[ , catCovarIndex[i]])[1]
      levelList[[i]] <- levels(tempDat[ , catCovarIndex[i]])
      #totalLevels <- sum(sapply(levelList, FUN = length))
    }
    catRefs <- unlist(refList)
  }else{
    catCovar <- FALSE
    catRefs <- NULL
    #totalLevels <- 0
  }



  #Create model formula
  #First consider the baseline covariates
  if(catCovar + contCovar > 0){
    baseString <- paste(paste0(" + Protein:", colnames(tempDat)[c(contIndex, catCovarIndex)]), collapse = "")
  }else(
    baseString <- ""
  )

  degVec <- as.character(1:timeDegree)
  degVec[1] <- ""

  if(randEffect == 0){
    fmla <- as.formula(paste("FC ~ 0 + Protein + Protein:Category", baseString, "+",
                             paste0("Protein:", paste0("Time", degVec), collapse = " + "),
                             "+", paste0("Protein:Category:", paste0("Time", degVec), collapse = " + ")))
  }
  if(randEffect == 1){
    fmla <- as.formula(paste("FC ~ 0 + Protein + Protein:Category", baseString, "+",
                             paste0("Protein:", paste0("Time", degVec), collapse = " + "),
                             "+", paste0("Protein:Category:", paste0("Time", degVec), collapse = " + ")))
    fmla <- paste(fmla, " + (1 | ", colnames(tempDat)[randIndex], ")")
  }
  if(randEffect == 2){
    degVec <- degVec[1] #force only linear slopes if random slopes are being fit to each ID
    fmla <- as.formula(paste("FC ~ 0 + Protein + Protein:Category", baseString, "+",
                             paste0("Protein:", paste0("Time", degVec), collapse = " + "),
                             "+", paste0("Protein:Category:", paste0("Time", degVec), collapse = " + ")))
    fmla <- paste(fmla, " + (1 + Time | ", colnames(tempDat)[randIndex], ")")
  }

  #Set reference category
  if(!is.null(refCat)){
    tempDat$Category <- factor(tempDat$Category)
    tempDat$Category <- relevel(tempDat$Category, ref = refCat)
  }else{
    tempDat$Category <- factor(tempDat$Category)
  }

  uCats <- levels(tempDat$Category)

  #create results matrix
  #rows = nProt.  cols = fulltimes*fullCats + 2*fullCats + levels(catCovar) + length(contIndex)
  #nPreds <- length(fullTimes) * length(fullCats)
  #nTests <- 2 * length(fullCats)
  #nCovars <- 2 * (length(contIndex) + totalLevels - length(catCovarIndex)) #estimate and pVal for each
  #deprecated

  tempNames <- paste0(rep(paste0("category:", fullCats), each = length(fullTimes) + 2),
                            c(paste0(":Time", fullTimes),"Pval-Time", "Pval-Category"))

  if(catCovar){
    baseCatNames <- unlist(lapply(levelList, function(x) x[-1]))
    catColNames <- paste0(c("Est_", "Pval_", "LL_", "UL_"), rep(baseCatNames, each = 4))
  }else{
    catColNames <- NULL
  }

  if(contCovar){
    baseContNames <- substring(colnames(tempDat)[contIndex], 22)
    contColNames <- paste0(c("Est_", "Pval_", "LL_", "UL_"), rep(baseContNames, each = 4))
  }else{
    contColNames <- NULL
  }

  tempNames <- c(tempNames, catColNames, contColNames)



  #Fit the model
  if(randEffect == 0){
    if(useW){
      fullMod <- lm(fmla, weights = w_, data = tempDat)
    }else{
      fullMod <- lm(fmla, data = tempDat)
    }
  }else{
    if(useW){
      fullMod <- lmer(fmla, weights = w_, data = tempDat)
    }else{
      fullMod <- lmer(fmla, data = tempDat)
    }
  }
  modSumm <- summary(fullMod)


  #The model fit may contain aliased variables
  #Prior to hypothesis testing, we need to know what Categories
  #are present for each protein
  obsCats <- by(tempDat, tempDat$Protein, FUN = function(x) which(fullCats %in% unique(x$Category)))
  refPos <- which(fullCats == refCat)


  ###############Now get predictions and p-values for each protein##############
  resMat <- matrix(NA, nrow = nProt, ncol = length(tempNames))
  colnames(resMat) <- tempNames

  for(index in 1:nProt){
    #First make sure the reference was observed in this protein
    if(!(refPos %in% obsCats[[index]])){
      next
    }else{
      refI <- which(obsCats[[index]] == refPos)
      dropRef <- obsCats[[index]][-refI]
    }

    #Implement F test for time effect
    #Create strings that define the hypothesis tests

    timeTests <- list()
    catTests <- list()


    #Test for an overall time effect in the baseline condition
    timeStr <-  paste0("Protein", uProt[index], paste0(":Time", degVec), " = 0")
    timeTests[[1]] <-  lht(fullMod, timeStr, singular.ok= T)$`Pr(>F)`[2]

    if(length(dropRef) > 0){

      for(t_ in 1:length(dropRef)){
        #Test for any difference between refCat and category t_
        catStr <-  paste0("Protein", uProt[index], ":Category", fullCats[dropRef[t_]], c("", paste0(":Time", degVec)), " = 0")
        catTests[[t_]] <- try(lht(fullMod, catStr, singular.ok = T)$`Pr(>F)`[2])

        #Test for an overall time effect in condition t_
        timeStr <- paste0("Protein", uProt[index], paste0(":Time", degVec), " + ",
                          "Protein", uProt[index], ":Category", fullCats[dropRef[t_]], paste0(":Time", degVec)," = 0")
        timeTests[[t_ + 1]] <- try(lht(fullMod, timeStr, singular.ok = T)$`Pr(>F)`[2])
      }#end condition loop

    }#end "if" more than one condition

    #if any tests failed something weird happened with aliasing
    catError <- unlist(lapply(catTests, function(x) attr(x, "class") == "try-error"))
    timeError <- unlist(lapply(timeTests, function(x) attr(x, "class") == "try-error"))
    if(length(catError) + length(timeError) > 0){next}

    #Figure out which times to predict
    #create list of times within observed ranges
    times <- list()
    protDat <- tempDat[which(tempDat$Protein == uProt[index]), ]
    for(c_ in 1:length(obsCats[[index]])){
      catDat <- protDat[which(protDat$Category == fullCats[obsCats[[index]][c_]]), ]
      #Check for missing values
      if(sum(is.na(catDat$FC)) > 0){
        catDat <- catDat[-which(is.na(catDat$FC)), ] #Remove missing values
      }

      catTime <- unique(catDat[ , "Time"])
      #catTime is the actual times observed for this protein within this category
      maxT <- max(catTime)
      minT <- min(catTime)
      boolT <- (fullTimes >= minT) & (fullTimes <= maxT)
      smalltime <- fullTimes[boolT]
      times[[c_]] <- smalltime
    }

    #Now extract and store the predicted values
    newDfs <- list()

    newDfs <- lapply(1:length(obsCats[[index]]), function(x)
        makePredDat(prot = uProt[index], timeVec = times[[x]], category = fullCats[obsCats[[index]][x]],
                    header = colnames(tempDat), timeDegree = timeDegree, catRefs = catRefs))


    catPreds <- lapply(newDfs, function(x) suppressWarnings(predict(fullMod, x)))

    tempPred <- matrix(c(catPreds[[refI]][match(fullTimes, times[[refI]])], timeTests[[1]]), nrow = 1)
    startPoint <- (refPos - 1) * (length(fullTimes) + 2) + 1
    resMat[index, startPoint:(startPoint + ncol(tempPred) - 1)] <- tempPred

    if(length(dropRef) > 0){
      for(t_ in 1:length(dropRef)){
        #need to map from dropRef to obsRef
        obsPos <- which(obsCats[[index]] == dropRef[t_])
        tempPred <- matrix(c(catPreds[[obsPos]][match(fullTimes, times[[obsPos]])],
                             timeTests[[t_ + 1]], catTests[[t_]]), nrow = 1)
        startPoint <- (dropRef[t_] - 1) * (length(fullTimes) + 2) + 1
        resMat[index, startPoint:(startPoint + ncol(tempPred) - 1)] <- tempPred

      }
    }


    #Now add point estimates, standard errors and pVals for each baseline covariate
    if(contCovar + catCovar > 0){
      #categorical comes before continuous in results table
      if(catCovar){
        startPoint <- (length(fullTimes) + 2) * length(fullCats) + 1
        for(k in 1:length(catCovarIndex)){ #loop through categorical variables
          catLevel <- levels(tempDat[ , catCovarIndex[k]])
          for(l in 1:(length(catLevel) - 1)){ #loop through levels of each variable
            levelName <- catLevel[l + 1]
            paramStr <- paste0("Protein", uProt[index], ":", colnames(tempDat)[catCovarIndex[k]],
                               levelName)
            #Make sure this level exists in the model
            coefIndex <- which(rownames(modSumm$coefficients) == paramStr)
            if(length(coefIndex) == 0){

              startPoint <- startPoint + 4
              next
                                      }
            tempRes <- matrix(c(modSumm$coefficients[coefIndex, c(1,4)],
                                confint(fullMod, paramStr)), nrow = 1)
            resMat[index, startPoint:(startPoint + 3)] <- tempRes
            startPoint <- startPoint + 4
          }

        }
      }else{
        startPoint <- (length(fullTimes) + 2) * length(fullCats) + 1 #pass start point to continuous covariate
      }

      if(contCovar){
        for(k in 1:length(contIndex)){
          paramStr <- paste0("Protein", uProt[index], ":", colnames(tempDat)[contIndex[k]])
          tempRes <- matrix(c(modSumm$coefficients[paramStr, c(1,4)],
                              confint(fullMod, paramStr)), nrow = 1)
          resMat[index, startPoint:(startPoint + 3)] <- tempRes
          startPoint <- startPoint + 4
        }
      }



      #baseDf <- do.call(cbind, baseRes)
      #pRes[[index]] <- cbind(pRes[[index]], baseDf)

    }#End addition of baseline covariate results


  }#end Gene loop

  #add identifiers
  if(groupByGene){
    resDf <- data.frame(Gene = tempDat$Gene[match(unique(tempDat$Protein), tempDat$Protein)],
                        Protein = savedProts[match(unique(tempDat$Protein), tempDat$Protein)], resMat)
  }else{
    resDf <- data.frame(Gene = tempDat$Gene[match(unique(tempDat$Protein), tempDat$Protein)],
                        Protein = unique(savedProts), resMat)
  }

  #Now add columns for random effects - ill conceived.  Return to this later
#  if(randEffect > 0){
 #   resDf <- cbind(resDf, ranef(fullMod))
  #}

  #Now add columns with q values
  #first find the p-values                      
  #pIndex <- grep("Pval", colnames(resDf))
  #Qvals <- as.data.frame(lapply(resDf[ , pIndex], function(x) p.adjust(x, method = "fdr")))
  #colnames(Qvals) <- gsub("Pval", "Qval", colnames(resDf)[pIndex])
  #finalDf <- data.frame(resDf, Qvals)

  #newIndex <- c(seq_along(resDf), pIndex + 0.5)
  #finalDf <- finalDf[ ,order(newIndex)]

  #finalDf
                       
  resDf
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
    testl <- lht(fullMod, paste(paste0("Protein", uProt[index], "_", paste0("Time", degVec), " = 0")))
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
