#File for functions that create plots

#' Caterpillar plot
#'
#' This function takes a list which contains a Stanfit object and summary
#' tables.  The output is a caterpillar plot of parameters with credible
#' intervals that do not contain 0.  Red line segments denote a 80% credible
#' interval, while the black tails show 95%.
#'
#' @export
#' @param results The results list generated from the function
#'   \code{\link{compCall}}
#' @param ptm A variable that determines what type of fold changes will be
#'   plotted.  The default, ptm = 0, plots relative protein abundance.
#'   Setting ptm = x, for any x > 0, will plot the PTM with ID x in the
#'   original dataframe. This will not work if x > 9.
#' @param allPars A boolean variable indicating whether or not all parameters
#'   should be plotted.  The default is false, resulting in a plot for only
#'   parameters with 95% credible intervals that do not contain zero.
#' @param byCond A boolean parameter which determines if separate plots should
#'   be made for each condition.  The default is FALSE which results in the
#'   creation of a single plot.
#'
caterpillar <- function(results, ptm = 0, allPars = FALSE,
                        byCond = FALSE){
  if(byCond == FALSE){
    if(ptm > 0){
      ptmName <- as.character(results[[2]]$ptmName)
      ptmType <- substring(ptmName, nchar(ptmName))
      ptmIndex <- which(ptmType == ptm)
      parSumm <- rstan::summary(results[[3]], pars = "alpha")$summary
      parName <- "alpha"
    }else{
      if(results[[3]]@par_dims$avgCond == 0){
        parSumm <- rstan::summary(results[[3]], pars = "beta")$summary
        parName <- "beta"
      }else{
        parSumm <- rstan::summary(results[[3]], pars = "avgCond")$summary
        parName <- "avgCond"
      }
    }

  if(allPars == FALSE){
    sigParIndex <- which(parSumm[ , "2.5%"] > 0 | parSumm[ , "97.5%"] < 0)
    if(ptm > 0){
      sigParIndex <- intersect(sigParIndex, ptmIndex)
    }
  }else{
    sigParIndex <- 1:dim(parSumm)[1]
    if(ptm > 0){
      sigParIndex <- intersect(sigParIndex, ptmIndex)
    }
  }
  orderedIndex <- sigParIndex[order(parSumm[sigParIndex, "50%"])]

  parStr <- paste(parName, "[", orderedIndex, "]", sep = "")

  rstan::plot(results[[3]], pars = parStr, mapping = ggplot2::theme) +
    ggplot2::scale_y_discrete(labels = NULL)
  }else{ # stratify by condition

    if(ptm > 0){
      ptmName <- as.character(results[[2]]$ptmName)
      ptmType <- substring(ptmName, nchar(ptmName))
      ptmIndex <- which(ptmType == ptm)
      parSumm <- rstan::summary(results[[3]], pars = "alpha")$summary
      parName <- "alpha"
      condition <- getCond(results[[2]]$ptmName, ptm = TRUE)
    }else{
      if(results[[3]]@par_dims$avgCond == 0){
        parSumm <- rstan::summary(results[[3]], pars = "beta")$summary
        parName <- "beta"
        condition <- getCond(results[[1]]$name, ptm = F)
      }else{
        parSumm <- rstan::summary(results[[3]], pars = "avgCond")$summary
        parName <- "avgCond"
        condition <- getCond(results[[1]]$name, ptm = F)
      }
    }

    for(i in 1:length(unique(condition))){
    condIndex <- which(condition == unique(condition)[i])

    if(allPars == FALSE){
      sigParIndex <- which(parSumm[ , "2.5%"] > 0 | parSumm[ , "97.5%"] < 0)
    }else{
      sigParIndex <- 1:dim(parSumm)[1]
    }
    #intersect the indices
    sigParIndex <- intersect(condIndex, sigParIndex)
    if(ptm > 0){
      sigParIndex <- intersect(sigParIndex, ptmIndex)
    }
    if(length(sigParIndex) == 0){next}

    orderedIndex <- sigParIndex[order(parSumm[sigParIndex, "50%"])]

    parStr <- paste(parName, "[", orderedIndex, "]", sep = "")

    cPlot <- rstan::plot(results[[3]], pars = parStr, mapping =
                           ggplot2::theme) +
      ggplot2::scale_y_discrete(labels = NULL) +
      ggplot2::ggtitle(paste("Condition", unique(condition)[i]))
    print(cPlot)
    } # end for loop

  } # end by condition = T

}#end caterpillar

#' Precision plot
#'
#' This function takes a results summary from the \code{\link{compCall}}
#' function.
#' The output is a plot of showing posterior means on the
#' x-axis and precision on the y-axis, where precision defined as the inverse
#' of the posterior variance.
#'
#' @export
#' @param summary A results table generated from the function
#'   \code{\link{compCall}}, which will be found in either the first (for
#'   protein summaries) or second (for ptm's) component of the result list
#'   object.
#' @param byCond A boolean parameter which determines if separate plots should
#'   be made for each condition.  The default is FALSE which results in the
#'   creation of a single plot.
#'
precisionPlot <- function(RES, ptm = 0, byCond = FALSE){
  colorVals <- c("(0,0.05]" = "#fb0000", "(0.05,0.25]" = "#dd1c77",
                                                 "(0.25,0.5]" = "#c994c7",
                                                 "(0.5,1]" = "black")
  labelScheme <- ggplot2::labs(y = "Log10 Precision", x = "Posterior mean of log2 fold-change")
  
  if(byCond){
    if(ptm == 0){
    condition <- getCond(RES[[1]]$name)
    sigFac <- cut(RES[[1]]$P_null, c(0,.05,.25,.5,1))
    newDf <- data.frame(RES[[1]], condition, "P_Null" = sigFac )
    ggplot2::ggplot(newDf, ggplot2::aes(x = mean, y = log10(1/var),
                                        shape = P_Null)) +
      ggplot2::geom_point(ggplot2::aes(color = P_Null)) + 
      ggplot2::scale_colour_manual(values = 
                                     c("(0,0.05]" = "#fb0000", 
                                       "(0.05,0.25]" = "#dd1c77",
                                       "(0.25,0.5]" = "#c994c7",
                                       "(0.5,1]" = "black")) + 
      labelScheme + ggplot2::facet_wrap( ~ condition, ncol = 5)
    }else{
      ptmType <- substring(RES[[2]]$ptmName, nchar(RES[[2]]$ptmName))
      ptmIndex <- which(ptmType == ptm)
      condition <- getCond(RES[[2]]$ptmName, TRUE)
      sigFac <- cut(RES[[2]]$P_null, c(0,.05,.25,.5,1))
      newDf <- data.frame(RES[[2]][ptmIndex, ], 
                          condition = condition[ptmIndex],
                          "P_Null" = sigFac[ptmIndex] )

      ggplot2::ggplot(newDf, ggplot2::aes(x = mean, y = log10(1/var))) +
        ggplot2::geom_point(ggplot2::aes(color = P_Null, shape = P_Null)) + 
        ggplot2::scale_colour_manual(values = 
                                       c("(0,0.05]" = "#fb0000", 
                                         "(0.05,0.25]" = "#dd1c77",
                                         "(0.25,0.5]" = "#c994c7",
                                         "(0.5,1]" = "black")) + 
        labelScheme + ggplot2::facet_wrap( ~ condition, ncol = 5)
    }
  }else{
    if(ptm == 0){
      sigFac <- cut(RES[[1]]$P_null, c(0,.05,.25,.5,1))
      newDf <- data.frame(RES[[1]], "P_Null" = sigFac )
      ggplot2::ggplot(newDf, ggplot2::aes(x = mean, y = log10(1/var))) +
      ggplot2::geom_point(ggplot2::aes(color = P_Null, shape = P_Null)) + 
      ggplot2::scale_colour_manual(values = 
                                     c("(0,0.05]" = "#fb0000", 
                                       "(0.05,0.25]" = "#dd1c77",
                                       "(0.25,0.5]" = "#c994c7",
                                       "(0.5,1]" = "black")) + labelScheme
    }else{
      ptmType <- substring(RES[[2]]$ptmName, nchar(RES[[2]]$ptmName))
      ptmIndex <- which(ptmType == ptm)
      sigFac <- cut(RES[[2]]$P_null, c(0,.05,.25,.5,1))
      newDf <- data.frame(RES[[2]], "P_Null" = sigFac )
      newDf <- newDf[ptmIndex, ]
      ggplot2::ggplot(newDf, ggplot2::aes(x = mean, y = log10(1/var))) +
        ggplot2::geom_point(ggplot2::aes(color = P_Null, shape = P_Null)) + labelScheme + ggplot2::scale_color_manual(values = 
                                        c("(0,0.05]" = "#fb0000", 
                                          "(0.05,0.25]" = "#dd1c77",
                                          "(0.25,0.5]" = "#c994c7",
                                          "(0.5,1]" = "black"))
    }
  }
}

#'
#' Variance plot
#'
#' This function shows the distribution of each parameter corresponding to an
#' experimental error.  Different variance components are modeled for each
#' tag/plex combination (along with different components for each PTM). If one
#' parameter is greatly different from the others it could be indicative of
#' problems in the experimental workflow.
#'
#' @export
#' @param model A results list generated from the function
#' \code{\link{compCall}}.
#'
checkVariance <- function(results){
  npars <- length(results[[4]]$varNames)
  ptmType <- substring(results[[4]]$varNames, nchar(results[[4]]$varNames))
  orderPars <- order(ptmType)
  parString <- paste("sigma[", orderPars, "]", sep = "")
  rstan::stan_plot(results[[3]], pars = parString) +
    ggplot2::scale_y_continuous(breaks = npars:1, 
                                labels = results[[4]]$varNames[orderPars])
}

