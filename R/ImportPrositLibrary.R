#' Loads a predicted spectral library and extracts the relevant information such as retention time, relative intensities of fragment ions etc.
#' @param zipfile absoute path to the zip file containing the spectral library file
#' @param filename name of the file to be extracted from the zip file
#' @export import.prosit.library

import.prosit.library <- function(zipfile,filename) {

  myPrositLib = read.csv(unz(zipfile,filename), quote="", stringsAsFactors=FALSE)
  print("Spectral library is loaded...")
  is2Z = myPrositLib[myPrositLib[,"PrecursorCharge"] == 2,]
  is3Z = myPrositLib[myPrositLib[,"PrecursorCharge"] == 3,]
  rm(myPrositLib)
  gc()

  startIndx2Z = which(!duplicated(is2Z[,"StrippedPeptide"]))
  endIndx2Z = c(startIndx2Z[2:length(startIndx2Z)]-1,nrow(is2Z))
  startIndx3Z = which(!duplicated(is3Z[,"StrippedPeptide"]))
  endIndx3Z = c(startIndx3Z[2:length(startIndx3Z)]-1,nrow(is3Z))
  print("Indices created for extraction of predicted spectra...")

  precursorCols = grep("StrippedPeptide|PrecursorMz|iRT|PrecursorCharge", colnames(is2Z)) # They are not sorted in the order of the pattern
  spectraCols = grep("RelativeIntensity|FragmentMz|FragmentType|FragmentNumber|FragmentCharge", colnames(is2Z)) # They are not sorted in the order of the pattern
  print("Extraction of relevant data columns...")
  predList = list("+2" = list(Precursor = cbind(startIndx2Z,endIndx2Z,is2Z[startIndx2Z,precursorCols]),Spectra = is2Z[,spectraCols]),
                  "+3" = list(Precursor = cbind(startIndx3Z,endIndx3Z,is3Z[startIndx3Z,precursorCols]),Spectra = is3Z[,spectraCols]))
  # Consider rewriting the code. It is unncecessary to create predList. It is better to create predOut directly instead!!!
  names(predList$`+2`$Precursor)[1:2] = c("startIndx","endIndx")
  names(predList$`+3`$Precursor)[1:2] = c("startIndx","endIndx")

  predList$`+3`$Precursor$startIndx = predList$`+3`$Precursor$startIndx + max(predList$`+2`$Precursor$endIndx)
  predList$`+3`$Precursor$endIndx = predList$`+3`$Precursor$endIndx + max(predList$`+2`$Precursor$endIndx)
  predOut = list(Precursor = rbind(predList$`+2`$Precursor,predList$`+3`$Precursor),
                 Spectra = rbind(predList$`+2`$Spectra,predList$`+3`$Spectra))

  print("Removing redundant entries...")
  rm(is2Z,is3Z,startIndx2Z,startIndx3Z,endIndx2Z,endIndx3Z,predList)
  gc()
  predOut

}
