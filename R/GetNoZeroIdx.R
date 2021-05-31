#' Get row indices for transitions in predicted spectral libraries where all intensities are not zero
#' @param files input spectral library files in generic format (csv)
#' @export get.no.zero.idx

get.no.zero.idx <- function(files) {

  for(i in 1:length(files)) {

    print(files[i])
    tmpInt = read_csv(paste(files[i], sep = ""),  col_types = cols_only(RelativeIntensity = col_double()))
    if(i == 1) {
      intMat = tmpInt
      intLength = length(tmpInt)
    } else if (intLength == length(tmpInt)){
      intMat[,paste("RelativeIntensity_", as.character(i), sep = "")] = tmpInt
      rm(tmpInt)
      gc()
    } else {
      print("Dimensions do not match...terminating code execution!")
      stop()
    }
  }
  notZero = which(apply(intMat != 0, 1, any) == T) # Get indices where all relative intensities are not zero, this will reduce library size with more than 25%.
}




