#Code for small silent functions used in the data manipulation


#Function that converts intensities into log ratios
#mat should be a matrix of intensities
#missing values and values < 1 will be set to 1.
logRatio <- function(mat, ref_index){
  mat[is.na(mat)] <- 1
  mat[mat < 1] <- 1
  lMat <- log(mat)
  lrMat <- lMat[ , -ref_index, drop = F] - rowMeans(lMat[ , ref_index, drop = F])
  lrMat
}

#Function that makes a single groupID from the first 4 rows of df
makeHeader <- function(df, index){

  header <- paste("lr", colnames(df)[index], df[1, index],
                  df[2, index], df[3, index], df[4, index], sep = "qqqq")
  header
}

#function that takes a dataframe and returns a dataframe with unique ids
transformDat <- function(df, modelFit, plexNumber){
  n_ <- nrow(df)
  value_index <- grep("tag", colnames(df))
  ref_index <- which(df[1, value_index] == df[1, value_index[1]])
  normal_index <- setdiff(1:length(value_index), ref_index)

  lrMat <- logRatio(as.matrix(df[5:(n_), value_index]), ref_index)
  header <- makeHeader(df[ , value_index], normal_index)
  colnames(lrMat) <- header

  newDf <- data.frame(Protein = df[5:n_, ]$Protein, Peptide = df[5:n_, ]$Peptide,                             bioID = df[5:n_, ]$bioID, Covariate =
                        df[5:n_, ]$Covariate, lrMat, stringsAsFactors = F)

  melted <- reshape2::melt(newDf, id.vars = c("Protein", "Peptide", "bioID",
                                              "Covariate"),
                 value.name = "lr", variable.name = "header")



  separated <- stringr::str_split_fixed(as.character(melted$header), "qqqq",6)

  #Merge tag with tenplex info
  tag_plex <- paste(separated[ , 2], melted$bioID, plexNumber, sep = "_")

  #figure out if we are using a bioid from the header or from a column
  bioCol <- df$bioID[1]
  if(bioCol == 1){
    bioID <- paste(melted$Protein, melted$bioID, sep = "_")
  }else{bioID <- paste(melted$Protein, separated[, 4])}

  finalDat <- data.frame(techID = paste(melted$Protein, separated[, 3], sep = "_"),
                         tag_plex, bioID,
                         ptmID = paste(melted$Peptide, separated[ , 6], sep = "_"),
                  condID = paste(melted$Protein, separated[, 5], sep = "_"),
                  ptm = separated[ , 6], covariate = melted$Covariate, lr = melted$lr
                         )
  finalDat
}#end function transformDat

