#' Loads a predicted spectral library and extracts the relevant information such as retention time, relative intensities of fragment ions etc.
#' @param zipfile absoute path to the zip file containing the spectral library file
#' @param noZeroIntensity logical indicating if all transitions with zero intensities should be removed
#' @param output absolute path to processed library file, as a RData object
#' @export process.prosit.lib

process.prosit.lib <- function(zipfile, noZeroIntensity, output) {


  print("Importing Prosit spectral library...")
  precursorCols = c("StrippedPeptide", "LabeledPeptide", "PrecursorMz", "PrecursorCharge", "iRT")
  msmsCols = c("RelativeIntensity", "FragmentMz", "FragmentType", "FragmentNumber", "FragmentCharge")
  myPrositLib = read_csv(zipfile, col_types = cols_only(StrippedPeptide = col_character(),
                                                        LabeledPeptide = col_character(),
                                                        PrecursorMz = col_double(),
                                                        PrecursorCharge = col_integer(),
                                                        iRT = col_double(),
                                                        RelativeIntensity = col_double(),
                                                        FragmentMz = col_double(),
                                                        FragmentType = col_character(),
                                                        FragmentNumber = col_integer(),
                                                        FragmentCharge= col_integer()))
  gc(full = T)
  startIdx = which(!duplicated(str_c(myPrositLib$LabeledPeptide, "_", myPrositLib$PrecursorCharge)))
  nRows = nrow(myPrositLib)
  #endIdx = c(startIdx[2:length(startIdx)] - 1, nrow(myPrositLib))
  if(noZeroIntensity) {
    print("Removing ALL transitions with intensities equal to zero...")
    myPrositLib = myPrositLib[-which(myPrositLib$RelativeIntensity == 0),]
    row.names(myPrositLib) = 1:nrow(myPrositLib)
    gc(full = T)
  }
  else {
    print("Including transitions with intensities equal to zero...")
  }
  if(any(myPrositLib$StrippedPeptide == "MACMACMACMACMACMACMACMACMACMAC" & myPrositLib$PrecursorCharge == 1)) {
    print("Removes dummy sequence")
    myPrositLib = myPrositLib[-which(myPrositLib$StrippedPeptide == "MACMACMACMACMACMACMACMACMACMAC" & myPrositLib$PrecursorCharge == 1), ]
    gc(full = T)
  }
  print(paste("Saving processed library to: ",output,sep = ""))
  outLib = list(PrecursorData = myPrositLib[startIdx, precursorCols],
                MsmsData = myPrositLib[startIdx,msmsCols])

  rm(myPrositLib)
  gc(full = T)
  outLib[["PrecursorData"]]$startIdx = startIdx
  outLib[["PrecursorData"]]$endIdx = c(startIdx[2:length(startIdx)] - 1, nRows)
  save(outLib,
       file = output,
       version = NULL,
       ascii = FALSE,
       compress = T,
       safe = T)
  gc(full = T)
  print("Processing completed!")

}
