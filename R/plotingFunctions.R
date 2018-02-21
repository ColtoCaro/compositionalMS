#File for functions that create plots

#' Caterpillar plot
#'
#' This function takes a list which contains a Stanfit object and summary
#' tables.  The output is a caterpillar plot of parameters with credible
#' intervals that do not contain 0.  Red line segments denote a 80\% credible
#' interval, while the black tails show 95\%.
#'
#' @export
#' @param RES The results list generated from the
#'   function \code{\link{compBayes}}
#' @param byCond A boolean parameter which determines if separate plots should
#'   be made for each condition.  The default is FALSE which results in the
#'   creation of a single plot.
#' @param allPars A boolean variable indicating whether or not all parameters
#'   should be plotted.  The default is false, resulting in a plot for only
#'   parameters with 95\% credible intervals that do not contain zero.
#'
#'
catterPlot <- function(RES, byCond = FALSE, plotAll = FALSE, avgCond = FALSE){

  bigDf <- do.call(rbind, RES[4:length(RES)])
  if(avgCond){
    tempCol <- colnames(bigDf)[1:8]
    bigDf <- bigDf[ , -c(3:9)]
    colnames(bigDf) <- tempCol
  }

  #find signif
  sigIndex <- which(bigDf$LL95 > 0 | bigDf$UL95 < 0)
  if(plotAll){sigIndex <- 1:nrow(bigDf)}

  newDf <- bigDf[sigIndex, ]
  newDf <- newDf[order(newDf$mean), ]
  newDf$index <- 1:nrow(newDf)

  labelScheme <- ggplot2::labs(y = "Log2 fold-change", x = "")
  if(byCond){

    ggplot2::ggplot(newDf, ggplot2::aes(y = mean, x = index)) +
      labelScheme +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = LL95, ymax = UL95)) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = LL80, ymax = UL80,
                                            colour = "80% Credible Interval"))+
      ggplot2::geom_point() +
      ggplot2::facet_wrap( ~ condition, ncol = 5)

  }else{

    ggplot2::ggplot(newDf, ggplot2::aes(y = mean, x = index)) +
      labelScheme +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = LL95, ymax = UL95)) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = LL80, ymax = UL80,
                                            colour = "80% Credible Interval"))+
      ggplot2::geom_point()


  }
}#end catterPlot

#' Precision plot
#'
#' This function takes a results summary from the \code{\link{compBayes}}
#' function.
#' The output is a plot of showing posterior means on the
#' x-axis and precision on the y-axis, where precision ploted as the inverse
#' of the posterior coefficient of variation.
#'
#' @export
#' @param summary A results table generated from the function
#'   \code{\link{compBayes}}, which will be found in either the first (for
#'   protein summaries) or second (for ptm's) component of the result list
#'   object.
#' @param byCond A boolean parameter which determines if separate plots should
#'   be made for each condition.  The default is FALSE which results in the
#'   creation of a single plot.
#' @param nullSet An interval which will determine the color scheme of the
#'   plot.  The probability of being in the nullSet interval is approximated
#'   and points are colored according to categories of this probability.  If
#'   a researcher is only interested in fold-changes greater than 2 then this
#'   interval should be set to (-1, 1).
#'
precisionPlot <- function(RES, byCond = FALSE, nullSet = c(-1,1), avgCond = FALSE){

  bigDf <- do.call(rbind, RES[4:length(RES)])
  if(avgCond){
    tempCol <- colnames(bigDf)[1:8]
    bigDf <- bigDf[ , -c(3:9)]
    colnames(bigDf) <- tempCol
  }

  labelScheme <- ggplot2::labs(y = "1 / CV",
                               x = "Posterior mean of log2 fold-change")
  if(byCond){

    pNull <- pnorm(nullSet[2], bigDf$mean, sqrt(bigDf$var)) -
      pnorm(nullSet[1], bigDf$mean, sqrt(bigDf$var))
    sigFac <- cut(pNull, c(0,.05,.25,.5,1), include.lowest = TRUE)

    cv <- sqrt(bigDf$var) / abs(bigDf$mean)

    newDf <- data.frame(bigDf, cv, pNull, "P_Null" = sigFac)
    ggplot2::ggplot(newDf, ggplot2::aes(x = mean, y = 1 / cv)) +
      ggplot2::geom_point(ggplot2::aes(color = P_Null))  +
      ggplot2::scale_colour_manual(values =
                                     c("[0,0.05]" = "#fb0000",
                                       "(0.05,0.25]" = "#dd1c77",
                                       "(0.25,0.5]" = "#c994c7",
                                       "(0.5,1]" = "black")) +
      labelScheme + ggplot2::facet_wrap( ~ condition, ncol = 5)

  }else{

    pNull <- pnorm(nullSet[2], bigDf$mean, sqrt(bigDf$var)) -
      pnorm(nullSet[1], bigDf$mean, sqrt(bigDf$var))
    sigFac <- cut(pNull, c(0,.05,.25,.5,1), include.lowest = TRUE)

    cv <- sqrt(bigDf$var) / abs(bigDf$mean)

    newDf <- data.frame(bigDf, cv, pNull, "P_Null" = sigFac)
    ggplot2::ggplot(newDf, ggplot2::aes(x = mean, y = 1 / cv)) +
      ggplot2::geom_point(ggplot2::aes(color = P_Null))  +
      ggplot2::scale_colour_manual(values =
                                     c("[0,0.05]" = "#fb0000",
                                       "(0.05,0.25]" = "#dd1c77",
                                       "(0.25,0.5]" = "#c994c7",
                                       "(0.5,1]" = "black")) +
      labelScheme


  }
}# end precisionPlot


#'
#' Variance plot
#'
#' This function shows the distribution of each parameter corresponding to an
#' experimental error.  Different variance components are modeled for each
#' tag/plex combination (along with different components for each PTM). If one
#' parameter is greatly different from the others it could be indicative of
#' problems in the experimental workflow.
#'
#' @param model A results list generated from the function
#' \code{\link{compBayes}}.
#'
checkVariance <- function(results){
  npars <- length(results[[4]]$varNames)
  #check to see if the model uses redundancies
  redun <- (length(grep("R",results[[4]]$varNames[1])) > 0)
  if(redun){
    rPos <- regexpr("R", results[[4]]$varNames)
    redundancy <- substring(results[[4]]$varNames, rPos + 1)
    ptmPos <- rPos - 1
  }else{ptmPos <- nchar(results[[4]]$varNames)}


  if(redun == F){
    ptmType <- substring(results[[4]]$varNames, ptmPos)
    orderPars <- order(ptmType)
    parString <- paste("sigma[", orderPars, "]", sep = "")
  rstan::stan_plot(results[[3]], pars = parString) +
    ggplot2::scale_y_continuous(breaks = npars:1,
                                labels = results[[4]]$varNames[orderPars])
  }else{
    rLevels <- levels(factor(redundancy))
    nPlots <- length(rLevels)
    for(i in 1:nPlots){
      whichVar <- which(redundancy == rLevels[i])
      parString <- paste("sigma[", whichVar, "]", sep = "")
      npars <- length(whichVar)
      pic <- rstan::plot(results[[3]], pars = parString) +
        ggplot2::scale_y_continuous(breaks = npars:1,
          labels = results[[4]]$varNames[whichVar])
      print(pic)
    }
  }


}#end checkVariance

