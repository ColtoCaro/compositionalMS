#Code for small silent functions used in the data manipulation


#Function that converts intensities into log ratios
#mat should be a matrix of intensities
#missing values and values < 1 will be set to 1.
logRatio <- function(mat, ref_index){
  mat[is.na(mat)] <- 1
  mat[mat < 1] <- 1
  lMat <- log2(mat)
  lrMat <- lMat[ , -ref_index, drop = F] - rowMeans(lMat[ , ref_index, drop = F])
  lrMat
}

#Function that makes a single groupID from the first 4 rows of df
makeHeader <- function(df, index){

  header <- paste("lr", colnames(df)[index], df[1, index],
                  df[2, index], df[3, index], sep = "qqqq")
  header
}

#function that takes a dataframe and returns a dataframe with unique ids
transformDat <- function(df, modelFit, plexNumber){
  n_ <- nrow(df)
  value_index <- grep("tag", colnames(df))
  condBio <- paste(df[1, value_index], df[2, value_index])
  ref_index <- which(condBio == condBio[1])
  normal_index <- setdiff(1:length(value_index), ref_index)

  lrMat <- logRatio(as.matrix(df[4:(n_), value_index]), ref_index)
  header <- makeHeader(df[ , value_index], normal_index)
  colnames(lrMat) <- header

  newDf <- data.frame(Protein = df[4:n_, ]$Protein, Peptide = df[4:n_, ]$Peptide,                             bioID = df[4:n_, ]$bioID, Covariate =
                        df[4:n_, ]$Covariate, lrMat, stringsAsFactors = F)

  melted <- reshape2::melt(newDf, id.vars = c("Protein", "Peptide", "bioID",
                                              "Covariate"),
                 value.name = "lr", variable.name = "header")



  separated <- stringr::str_split_fixed(as.character(melted$header), "qqqq",5)

  #Merge tag with tenplex info
  tag_plex <- paste(separated[ , 2], melted$bioID, plexNumber, sep = "_")

  #figure out if we are using a bioid from the header or from a column
  bioCol <- df$bioID[1]
  if(bioCol == 1){
    bioID <- paste(melted$Protein, separated[ , 3], melted$bioID,  sep = "_")
  }else{bioID <- paste(melted$Protein, separated[ , 3], separated[, 4])}

  finalDat <- data.frame(condID = paste(melted$Protein, separated[, 3],
                                        sep = "_"), bioID,
                         ptmID = paste(melted$Protein, separated[, 3],
                                       melted$Peptide, separated[ , 5],
                                       sep = "_"),
                  ptm = separated[ , 5], tag_plex, covariate = melted$Covariate,
                  lr = melted$lr, stringsAsFactors = F)
  finalDat <- finalDat[order(finalDat$tag_plex, finalDat$condID, finalDat$bioID, finalDat$ptm, finalDat$ptmID), ]

  #normalize the non-ptm data by tag
  normed <- unlist(by(melted[finalDat$ptm == 0, ]$lr,
                 finalDat[finalDat$ptm == 0, ]$tag_plex, function(x)
                   x - mean(x)))
  melted$lr[finalDat$ptm == 0] <- normed
  finalDat$lr <- melted$lr

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


