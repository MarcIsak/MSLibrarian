#' Extracts the unique precursor information from an imported calibration library
#' @param calibrationLib input Calibration library object
#' @export get.lib.precursors2

get.lib.precursors2 <- function(calibrationLib) {

  preCols = c("PrecursorMz", "PrecursorCharge", "NormalizedRetentionTime",
              "PeptideSequence", "ModifiedPeptideSequence", "ProteinId", "precursors")
  print("Extracting precursor data...")
  rawMat = calibrationLib@RawLib
  rawMat$precursors = str_c(rawMat$ModifiedPeptideSequence, "_", rawMat$PrecursorCharge)
  startIdx = which(!duplicated(rawMat$precursors))
  endIdx = c(startIdx[2:length(startIdx)]-1, nrow(rawMat))
  print("Adding precursor data to the Calibration Library object...")
  calibrationLib@PrecursorData$Original = rawMat[startIdx,preCols]
  calibrationLib@PrecursorData$Original$startIdx = as.integer(startIdx)
  calibrationLib@PrecursorData$Original$endIdx = as.integer(endIdx)
  calibrationLib@PrecursorData$Original$msmsIdx = seq(1,nrow(calibrationLib@PrecursorData$Original))
  calibrationLib@PrecursorData$Original$precursors = NULL
  calibrationLib

}
