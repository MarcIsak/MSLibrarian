#' A filter which sets the maximum number of of transitions that a library entry can have
#' @param msLib an input MSLibrarian object
#' @param predLib character defining the predicted spectral library to filter c("PrositLib, "Ms2pipLib")
#' @param maxTransitions integer setting the maximum number of transitions of a library entry
#' @export max.transitions

max.transitions <- function(msLib, predLib, maxTransitions) {

  get.max.transitions <- function(idx, maxTransitions) {

    msmsData = slot(msLib@PredLib,predLib)$MsmsData[seq(min(as.numeric(idx)),max(as.numeric(idx))),]
    if(length(seq(min(as.numeric(idx)),max(as.numeric(idx)))) > maxTransitions) {
      msmsDataFilter = msmsData[order(msmsData$RelativeIntensity, decreasing = T)[1:maxTransitions],]
    } else {
      msmsData
    }
  }
  print(paste("Filtering library entries using a maximum of ", as.character(maxTransitions), " transitions!", sep = ""))
  libList = apply(slot(msLib@PredLib,predLib)$PrecursorData[,c("startIdx", "endIdx")], 1, get.max.transitions, maxTransitions)
  print("Creating new precursor data indices for MS/MS data")
  startIdx = c(1,cumsum(unlist(lapply(libList, nrow))) + 1)
  endIdx = startIdx[2:length(startIdx)]-1
  startIdx = startIdx[1:length(startIdx)-1]
  print("Adding filtered MS/MS data to slot")
  slot(msLib@PredLib,predLib)$MsmsData = do.call('rbind', libList)
  print("Setting new row name indices for MS/MS data")
  row.names(slot(msLib@PredLib,predLib)$MsmsData) = 1:nrow(slot(msLib@PredLib,predLib)$MsmsData)
  slot(msLib@PredLib,predLib)$PrecursorData[,c("startIdx", "endIdx")] = data.frame(startIdx = startIdx, endIdx = endIdx, stringsAsFactors = F)
  print("Processing completed!")
  msLib

}
