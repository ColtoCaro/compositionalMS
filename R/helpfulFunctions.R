#Code for small silent functions used in the data manipulation


#Function that returns the number of ampersands in column headers
checkQs <- function(dat){
  ampList <- lapply(dat, function(x) sum(grepl('qqqq', names(x), fixed=T)))
  n_amps <- sum(unlist(ampList))
  n_amps
}
checkNames <- function(dat){
  repList <- lapply(dat, function(x) sum(duplicated(colnames(x))))
  n_reps <- sum(unlist(repList))
  n_reps
}
#Function that converts intensities into log ratios
#mat should be a matrix of intensities
#missing values and values < 1 will be set to 1.
logRatio <- function(mat){
  mat[is.na(mat)] <- 1
  mat[mat < 1] <- 1
  lMat <- log(mat)
  lrMat <- lMat[ , -1] - lMat[ , 1]
  lrMat
}

#Function that makes a single groupID from the first 3 rows
makeHeader <- function(df){

  header <- paste("lr", colnames(df)[5:length(df)], df[1, 5:length(df)],
                  df[2, 5:length(df)], df[3, 5:length(df)], sep = "qqqq")
  header
}

#function that takes a dataframe and returns a dataframe with unique
#protein ids
transformDat <- function(df, modelFit){
  p_ <- length(df)
  n_ <- nrow(df)
  lrMat <- logRatio(as.matrix(df[4:(n_), 4:(p_)]))

  header <- makeHeader(df)
  colnames(lrMat) <- header
  newDf <- data.frame(Protein = df[4:n_, ]$Protein, Covariate =
                        df[4:n_, ]$Covariate, lrMat, stringsAsFactors = F)

  melted <- melt(newDf, id.vars = c("Protein", "Covariate"), value.name = "lr"
                 , variable.name = "header")

  separated <- stringr::str_split_fixed(as.character(melted$header), "qqqq",5)
  finalDat <- data.frame(techID = paste(melted$Protein, separated[, 2],
                                        separated[, 3]),
                  bioID = paste(melted$Protein, separated[, 2], separated[, 4]),
                  condID = paste(melted$Protein, separated[, 2], separated[, 5]),
                  covariate = melted$Covariate, lr = melted$lr
                         )

}#end function transformDat

