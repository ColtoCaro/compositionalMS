#Code for small silent functions used in the data manipulation

#Function that converts intensities into log ratios
#mat should be a matrix of intensities
logRatio <- function(mat){
  lMat <- log(mat)
  lrMat <- lMat[ , -1] - lMat[ , 1]
  lrMat
}

#Function that makes a single groupID from the first 3 rows
makeHeader <- function(df){

  header <- paste("lr", colnames(df)[5:length(df)], df[1, 5:length(df)],
                  df[2, 5:length(df)], df[3, 5:length(df)], sep = "&")
  header
}

#function that takes a dataframe and returns a dataframe with unique
#protein ids
transformDat <- function(df){
  p_ <- length(df)
  n_ <- nrow(df)
  lrMat <- logRatio(as.matrix(df[4:(n_), 4:(p_)]))

  header <- makeHeader(df)
  colnames(lrMat) <- header
  newDf <- data.frame(Protein = df[4:n_, ]$Protein, Covariate =
                        df[4:n_, ]$Covariate, lrMat, stringsAsFactors = F)

  melted <- melt(newDf, id.vars = c("Protein", "Covariate"), value.name = "lr"
                 , variable.name = "header")

  separated <- strsplit(as.character(melted$header), ".", fixed = T)

}#end function transformDat

