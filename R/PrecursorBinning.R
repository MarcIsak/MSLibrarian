#' Bin precursors based on number cutoff
#' @param ceMat input dot product matrix from predicted vs. experimental spectral matches.
#' @param cutoff minimum number of precursors to have a separate bin
#' @export precursor.binning

precursor.binning <- function(ceMat, cutoff) {

  set.range <- function(ranges, mat, z) {
    print(z)
    if(length(ranges) == 1) {
      ranges = c(ranges, ranges + 1)
      #print(which(between(mat$peptideLength, lower = min(ranges), upper = max(ranges)) & mat$charge == z))
      mat$lenRange[between(mat$peptideLength, lower = min(ranges), upper = max(ranges)) & mat$charge == z] = str_c(range(ranges), collapse = "_")
    } else if (length(ranges) > 1) {
      mat$lenRange[between(mat$peptideLength, lower = min(ranges), upper = max(ranges)) & mat$charge == z] = str_c(range(ranges), collapse = "_")
      #print(which(between(mat$peptideLength, lower = min(ranges), upper = max(ranges)) & mat$charge == z))
    }
    mat
  }
  print(str_c("Bin precursors based on a cutoff: ", cutoff))
  ceMat$lenRange = as.character(ceMat$peptideLength)

  for (i in as.numeric(unique(ceMat$charge))) {
    test = as.data.frame(table(ceMat$peptideLength[ceMat$charge == i & ceMat$ce == min(ceMat$ce)]), stringsAsFactors = F)
    colnames(test) = c("peptideLength", "nbr")
    test$peptideLength = as.numeric(test$peptideLength)
    #test$nbr = as.numeric(test$nbr)
    idx = which(test$nbr > cutoff)
    subPrec = test[idx,]

    ceMat = set.range(test$peptideLength[which(test$peptideLength < min(subPrec$peptideLength))], ceMat, i)
    ceMat = set.range(test$peptideLength[which(test$peptideLength > max(subPrec$peptideLength))], ceMat, i)
  }
  ceMat
}
