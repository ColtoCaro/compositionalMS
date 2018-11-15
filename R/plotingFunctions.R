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
#' @param plotAll A boolean variable indicating whether or not all parameters
#'   should be plotted.  The default is false, resulting in a plot for only
#'   parameters with 95\% credible intervals that do not contain zero.
#' @param avgCond A boolean parameter that determines whether or not plots
#'   will be based on heierarchical population level mean parameters or average
#'   within sample parameters.  The default is to use the population level
#'   parameter, however with very small sample sizes this plot may not be useful.
#'
#'
catterPlot <- function(RES, byCond = FALSE, plotAll = FALSE, avgCond = FALSE){
  if(avgCond){
    tempTab <- RES[[3]]
  }else{
    tempTab <- RES[[5]]
  }
  meltEst <- melt(tempTab[ , c(grep("Protein", colnames(tempTab)), grep("Est", colnames(tempTab)))],
                  id.vars = "Protein", variable.name = "condition", value.name = "mean")
  meltEst <- meltEst[order(meltEst$Protein), ]
  meltLL <- melt(tempTab[ , c(grep("Protein", colnames(tempTab)), grep("LL", colnames(tempTab)))],
                  id.vars = "Protein", variable.name = "condition", value.name = "LL95")
  meltLL <- meltLL[order(meltLL$Protein), ]
  meltUL <- melt(tempTab[ , c(grep("Protein", colnames(tempTab)), grep("UL", colnames(tempTab)))],
                  id.vars = "Protein", variable.name = "condition", value.name = "UL95")
  meltUL <- meltUL[order(meltUL$Protein), ]

  bigDf <- data.frame(meltEst, LL95 = meltLL$LL95, UL95 = meltUL$UL95)

  # bigDf <- do.call(rbind, RES[4:length(RES)])
  # if(avgCond){
  #   tempCol <- colnames(bigDf)[1:8]
  #   bigDf <- bigDf[ , -c(3:9)]
  #   colnames(bigDf) <- tempCol
  # }

  #find signif
  sigIndex <- which(bigDf$LL95 > 0 | bigDf$UL95 < 0)
  if(plotAll){sigIndex <- 1:nrow(bigDf)}
  if(length(sigIndex) == 0){stop("There are no significant changes.  Use plotAll = T to see results")}

  newDf <- bigDf[sigIndex, ]
  newDf <- newDf[order(newDf$mean), ]
  newDf$index <- 1:nrow(newDf)

  labelScheme <- ggplot2::labs(y = "Log2 fold-change", x = "")
  if(byCond){

    ggplot2::ggplot(newDf, ggplot2::aes(y = mean, x = index)) +
      labelScheme +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = LL95, ymax = UL95)) +
      #ggplot2::geom_pointrange(ggplot2::aes(ymin = LL80, ymax = UL80,
      #                                      colour = "80% Credible Interval"))+
      ggplot2::geom_point() +
      ggplot2::facet_wrap( ~ condition, ncol = 5)

  }else{

    ggplot2::ggplot(newDf, ggplot2::aes(y = mean, x = index)) +
      labelScheme +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = LL95, ymax = UL95)) +
      #ggplot2::geom_pointrange(ggplot2::aes(ymin = LL80, ymax = UL80,
      #                                      colour = "80% Credible Interval"))+
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
#' @param RES A results list generated from the function
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
#' @param avgCond A boolean parameter that determines whether or not plots
#'   will be based on heierarchical population level mean parameters or average
#'   within sample parameters.  The default is to use the population level
#'   parameter, however with very small sample sizes this plot may not be useful.
#'
precisionPlot <- function(RES, byCond = FALSE, nullSet = c(-1,1), avgCond = FALSE, ptm = FALSE){
  if(ptm == FALSE){
  if(avgCond){
    tempTab <- RES[[3]]
  }else{
    tempTab <- RES[[5]]
  }
  meltEst <- melt(tempTab[ , c(grep("Protein", colnames(tempTab)), grep("Est", colnames(tempTab)))],
                  id.vars = "Protein", variable.name = "condition", value.name = "mean")
  meltEst <- meltEst[order(meltEst$Protein), ]
  meltVar <- melt(tempTab[ , c(grep("Protein", colnames(tempTab)), grep("Var", colnames(tempTab)))],
                 id.vars = "Protein", variable.name = "condition", value.name = "var")
  meltVar <- meltVar[order(meltVar$Protein), ]

  }else{
    tempTab <- RES[[9]]

    meltEst <- melt(tempTab[ , c(grep("Peptide", colnames(tempTab)), grep("Est", colnames(tempTab)))],
                    id.vars = "Peptide", variable.name = "condition", value.name = "mean")
    meltEst <- meltEst[order(meltEst$Peptide), ]
    meltVar <- melt(tempTab[ , c(grep("Peptide", colnames(tempTab)), grep("Var", colnames(tempTab)))],
                    id.vars = "Peptide", variable.name = "condition", value.name = "var")
    meltVar <- meltVar[order(meltVar$Peptide), ]
  }

  bigDf <- data.frame(meltEst, var = meltVar$var)

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
                                       "(0.05,0.25]" = "#c994c7", #"#dd1c77",
                                       "(0.25,0.5]" =  "violetRed4",
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
                                       "(0.05,0.25]" = "#c994c7", #"#dd1c77",
                                       "(0.25,0.5]" = "violetRed4",
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


#' Proportion plot
#'
#' This function plots the proportion of ions belonging
#' to each biological replicate.  The first replicate
#' is typically a bridge channel and all repicates are
#' ordered and colored by experimental condition.
#'
#' @export
#' @param df A dataframe typically generated and stored in the second
#'   component of a results list from compBayes.
#' @param gene A string that determines what data will be plotted.
#' @param condKey A matrix with two columns.  The first column contains
#'   integers denoting the experimental condition for each biological
#'   replicate.  The second contains integers that represent each
#'   biological replicate.  This key should be taken from the
#'   header information required to run the compBayes function.
#'
propPlot <- function(df, gene = NULL, condKey, protName = NULL){
  if(is.null(gene) & is.null(protName)){
    stop("Please specify a gene or protein to plot")
  }

  estIndex <- grep("Est" , colnames(df))
  ulIndex <- grep("UL95" , colnames(df))
  llIndex <- grep("LL95", colnames(df))

  if(!is.null(protName)){
    geneIndex <- match(protName, df$Protein)
    plotName <- protName
  }else{
    geneIndex <- grep(gene, df$Gene)
    plotName <- gene
  }

  if(length(geneIndex) > 1){
     print("Warning: Gene name was not unique.  First discovered entry was plotted \n
           Specify protein name to see a different plot")
     geneIndex <- geneIndex[1]
  }

  channel <- colnames(df)[estIndex]
  bioRep <- as.integer(substring(channel, regexpr("D", channel) + 1))
  cond <- condKey[match(bioRep, condKey[ , 2]) , 1]

  newDat <- data.frame(
    Condition = cond,
    Channel = bioRep,
    Est = t(df[geneIndex, estIndex]),
    LL = t(df[geneIndex, llIndex]),
    UL = t(df[geneIndex, ulIndex])
  )
  names(newDat) <- c("Condition", "Channel", "Est", "LL", "UL")

  newDat <- newDat[order(newDat$Condition, newDat$Channel), ]
  newDat$Channel <- factor(newDat$Channel, levels = unique(newDat$Channel))
  newDat$Condition <- factor(newDat$Condition, levels = unique(newDat$Condition))

  titleStr <- paste0("Proportion Plot for ", plotName)
  labelScheme <- ggplot2::labs(y = "Estimated Proportion", x = "Biological Replicate")
  ggplot2::ggplot(newDat, ggplot2::aes(y = Est, x = Channel, colour = Condition)) +
    labelScheme +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = LL, ymax = UL)) +
    ggplot2::geom_point() +
    ggplot2::ggtitle(titleStr)


}


