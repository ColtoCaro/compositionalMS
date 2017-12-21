#Code for small silent functions used in the data manipulation


#Function that converts intensities into log ratios
#also returns the minimum log intensity for each pair
#mat should be a matrix of intensities
#missing values and values < 1 will be set to 1.
logRatio <- function(mat, ref_index){
  #mat[is.na(mat)] <- 1  normalization now done prior to this routine
  #mat[mat < 1] <- 1
  lMat <- log2(mat)
  denom <- rowMeans(lMat[ , ref_index, drop = F])
  numCols <- lMat[ , -ref_index, drop = F]

  p_ <- ncol(numCols)
  denomMat <- matrix(rep(denom, p_), ncol = p_)
  mintensity <- 2^(denomMat) + 2^(numCols) #for making pairwise ssn
  lrMat <- numCols - denom

  list(lrMat, mintensity)
}


#Function that makes a single groupID from the first 4 rows of df
makeHeader <- function(df, index){

  header <- paste("lr", colnames(df)[index], df[1, index],
                  df[2, index], df[3, index], sep = "qqqq")
  header
}

#function that takes a dataframe and returns a dataframe with unique ids
transformDat <- function(df, plexNumber, normalize){
  #convert factors to strings
  facIndex <- which(sapply(df, is.factor))
  df[facIndex] <- lapply(df[facIndex], as.character)

  #Zero out unused columns
  bioCol <- df$bioID[1]
  if(bioCol == 0){
    df$bioID[] <- 0
  }

  covarCol <- df$Covariate[1]
  if(covarCol == 0){
    df$Covariate[] <- 0
  }
  varCol <- df$varCat[1]
  if(varCol == 0){
    df$varCat[] <- 0
  }


  if(df[2, 1] == 0){
    df[2, ][] <- 0
  }
  if(df[3, 1] == 0){
    df[3, ][] <- 0
  }

  n_ <- nrow(df)

  jDat <- df[4:(n_), ]
  df[4:n_, ] <- jDat[order(jDat$bioID), ]

  value_index <- grep("tag", colnames(df))

  nMat <- df[4:(n_), value_index]
  #normalize the df by bioRep
  if(normalize == TRUE){
  normed <- by(as.matrix(df[4:(n_), value_index]), df$bioID[4:n_], cNorm)
  nMat <- as.matrix(do.call(cbind, normed))
  }

  condBio <- paste(df[1, value_index],
                   as.integer(df[2, 1])*df[2, value_index])
  ref_index <- which(condBio == condBio[1])
  normal_index <- setdiff(1:length(value_index), ref_index)

  lRes <- logRatio(nMat, ref_index)
  lrMat <- lRes[[1]]
  minTensities <- lRes[[2]]
  header <- makeHeader(df[ , value_index], normal_index)
  colnames(lrMat) <- header

  newDf1 <- data.frame(Protein = df[4:n_, ]$Protein,
                      Peptide = df[4:n_, ]$Peptide,                                                  bioID = df[4:n_, ]$bioID,
                      Covariate = df[4:n_, ]$Covariate,
                      varCat = df[4:n_, ]$varCat,
                      lrMat, stringsAsFactors = F)

  newDf2 <- data.frame(Protein = df[4:n_, ]$Protein,
                       Peptide = df[4:n_, ]$Peptide,                                                  bioID = df[4:n_, ]$bioID,
                       Covariate = df[4:n_, ]$Covariate,
                       varCat = df[4:n_, ]$varCat,
                       minTensities, stringsAsFactors = F)

  melted1 <- reshape2::melt(newDf1, id.vars = c("Protein", "Peptide", "bioID",
                                              "Covariate", "varCat"),
                 value.name = "lr", variable.name = "header")


  melted2 <- reshape2::melt(newDf2, id.vars = c("Protein", "Peptide", "bioID",
                                                "Covariate", "varCat"),
                            value.name = "pairMin", variable.name = "header")


  melted <- data.frame(melted1, pairMin = melted2$pairMin)

  separated <- stringr::str_split_fixed(as.character(melted$header), "qqqq",5)

  #Merge tag with tenplex info

  tag_plex <- paste(separated[ , 2], melted$bioID, plexNumber, sep = "_")

  #figure out if we are using a bioid from the header or from a column

  if(bioCol == 1){
    bioID <- paste(melted$Protein, separated[ , 3], melted$bioID,  sep = "_")
  }else{bioID <- paste(melted$Protein, separated[ , 3], separated[, 4])}

  finalDat <- data.frame(protein = melted$Protein,
                  condID = paste(melted$Protein, separated[, 3],
                                        sep = "_"), bioID,
                  ptmID = paste(melted$Protein, separated[, 3],
                                       melted$Peptide, separated[ , 5],
                                       sep = "_"),
                  ptm = separated[ , 5], tag_plex,
                  covariate = melted$Covariate,
                  varCat = melted$varCat,
                  pairMin = melted$pairMin,
                  lr = melted$lr, stringsAsFactors = F)
  finalDat <- finalDat[order(finalDat$tag_plex, finalDat$condID,
                             finalDat$bioID, finalDat$ptm, finalDat$ptmID), ]

  #do normalization at a previous step
  #normalize the non-ptm data by tag

  #normed <- unlist(by(melted[finalDat$ptm == 0, ]$lr,
   #              finalDat[finalDat$ptm == 0, ]$tag_plex, function(x)
    #               x - mean(x)))
  #melted$lr[finalDat$ptm == 0] <- normed
  #finalDat$lr <- melted$lr


  finalDat
}#end function transformDat

#function for extracting the condition number from labels
getCond <- function(strVec, ptm = FALSE){
  sPosition <- regexpr("_", strVec)
  if(ptm){
  subbed <- sub("_", "*", strVec)
  ePosition <- regexpr("_", subbed)
  condNumber <- as.integer(substring(strVec, sPosition +1, ePosition -1))
  }else{
    condNumber <- as.integer(substring(strVec, sPosition +1))
  }

  condNumber
}

#function to reverse strings
strReverse <- function(x){
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}

#closure function
cl <- function(vec){
  vec/sum(vec)
}

#centered log ratio
clr <- function(vec){
  gm <- exp(mean(log(vec)))
  log(vec/gm)
}

clrInv <- function(vec){
  cl(exp(vec))
}

#Function to normalize the data

cNorm <- function(mat, subIndex = NULL){
  #takes a compostion matrix and subracts out the mean
  mat[mat < 1] <- 1
  clrMat <- t(apply(mat, 1, clr))

  if(is.null(subIndex)){
    cMean <- clrInv(apply(clrMat, 2, mean))
  }else{
    cMean <- clrInv(apply(clrMat[subIndex, ], 2, mean))
  }

  normed <- t(apply(mat, 1, function(x) cl(x/cMean)))
  normed
}


#Function to prepare data for an ANOVA.

anovaPrep <- function(df){
  #Zero out unused column
  bioCol <- df$bioID[1]
  if(bioCol == 0){
    df$bioID[] <- 0
  }
  n_ <- nrow(df)

  jDat <- df[4:(n_), ]
  df[4:n_, ] <- jDat[order(jDat$bioID), ]

  value_index <- grep("tag", colnames(df))

header <- makeHeader(df[ , value_index], normal_index)
  colnames(lrMat) <- header

  newDf1 <- data.frame(Protein = df[4:n_, ]$Protein,
                       Peptide = df[4:n_, ]$Peptide,                                                  bioID = df[4:n_, ]$bioID,
                       Covariate = df[4:n_, ]$Covariate,
                       varCat = df[4:n_, ]$varCat,
                       lrMat, stringsAsFactors = F)

}

#function to make a matrix of relevant samples
makeMat <- function(vec, model){
  extracted <- lapply(vec, function(x)
    rstan::extract(model,
                   pars=paste("beta[", x, "]", sep = ""))$beta)
  extractMat <- do.call(cbind, extracted)
  extractMat
}

#function to apply alr inverse to a matrix
alrInv <- function(mat){

zeroMat <- cbind(rep(0,nrow(mat)),mat)
expMat <- 2^zeroMat
simplex <- t(apply(expMat, 1, function(x) x/sum(x)))
simpMeans <- apply(simplex, 2, mean)
simpVar <- apply(simplex, 2, var)
simpInt <- apply(simplex, 2, quantile, probs = c(.025, .975, .1, .9))

list(simpMeans, simpVar, simpInt)
}
