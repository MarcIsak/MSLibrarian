#' Loads a predicted spectral library and extracts the relevant information such as retention time, relative intensities of fragment ions etc.
#' @param zipfile absoute path to the zip file containing the spectral library file
#' @param noZeroIntensity logical indicating if all transitions with zero intensities should be removed
#' @param transitionFilter logical indicating if a transition filter should be used
#' @param transitionThreshold numeric value (0,1) setting the percentage of most observed transitions to include
#' @param output absolute path to processed library file, as a RData object
#' @export process.prosit.lib

process.prosit.lib <- function(zipfile, noZeroIntensity, transitionFilter, transitionThreshold, output) {

  prosit.list <- function(x) {
    myPrositLib[myPrositLib[,"PrecursorCharge"] == x,]
  }
  add.idx <- function(libList, offsetIdx) {

    cols = "StrippedPeptide|LabeledPeptide|PrecursorMz|iRT|PrecursorCharge"
    startIdx = which(!duplicated(str_c(libList$LabeledPeptide, "_", libList$PrecursorCharge)))
    # startIdx = which(libList$FragmentNumber == 1 &
    #                    libList$FragmentType == "y" &
    #                    libList$FragmentCharge == 1)
    outList = libList[startIdx, grep(cols, colnames(libList))]
    outList$startIdx = startIdx + offsetIdx[i]
    outList$endIdx = c(outList$startIdx[2:length(outList$startIdx)]-1, nrow(libList) + offsetIdx[i])
    i <<- i + 1
    outList
  }
  print("Importing Prosit spectral library...")
  myPrositLib = as.data.frame(read_csv(zipfile), stringsAsFactors = F)
  gc(full = T)
  if(noZeroIntensity) {
    print("Removing ALL transitions with intensities equal to zero...")
    myPrositLib = myPrositLib[-which(myPrositLib$RelativeIntensity == 0),]
    row.names(myPrositLib) = 1:nrow(myPrositLib)
    gc(full = T)
  }
  else {
    print("Including transitions with intensities equal to zero...")
  }
  if(transitionFilter) {

    idx = msLib@SpectrastLib@Unfiltered@PrecursorData$FilterLib[, c("startIdx", "endIdx")]

    annot = table(unlist(apply(idx, 1, get.annot)))
    sortAnnot = sort(annot, decreasing = T)/sum(annot)

    transitions = names(sortAnnot)[which(cumsum(sortAnnot) <= transitionThreshold)]

    prositAnnots = str_c(myPrositLib$FragmentType, myPrositLib$FragmentNumber, "^",
                         myPrositLib$FragmentCharge)
    selPrositAnnots = grepl(paste(str_replace(transitions, "\\^", "\\\\^"),
                                  collapse = "|"),
                            prositAnnots)

    myPrositLib = myPrositLib[selPrositAnnots,]
    print(paste("Keeping only",
                as.character(length(table(prositAnnots[selPrositAnnots]))),
                "/",
                as.character(length(table(prositAnnots))),
                "transition types in the predicted library"))

  } else {
    prositAnnots = str_c(myPrositLib$FragmentType, myPrositLib$FragmentNumber, "^",
                         myPrositLib$FragmentCharge)
    print(paste("Keeping all", as.character(length(table(prositAnnots))),"transition types in the predicted library"))
  }
  if(any(myPrositLib$StrippedPeptide == "MACMACMACMACMACMACMACMACMACMAC" & myPrositLib$PrecursorCharge == 1)) {
    print("Removes dummy sequence")
    myPrositLib = myPrositLib[-which(myPrositLib$StrippedPeptide == "MACMACMACMACMACMACMACMACMACMAC" & myPrositLib$PrecursorCharge == 1), ]
    gc(full = T)
  }
  i <<- 1
  print("Sorting spectral library by precursor charges...")
  libList = lapply(unique(myPrositLib$PrecursorCharge),prosit.list)
  rm(myPrositLib)
  gc(full = T)
  print("Indexing of spectral library...")
  print("Removal of redundant entries...")
  i <<- 1
  offsetIdx = c(0, unlist(lapply(libList, nrow)))
  msmsCols = "RelativeIntensity|FragmentMz|FragmentType|FragmentNumber|FragmentCharge"

  # Remove all in this section later --

  #save(libList, file = output, version = NULL, ascii = FALSE, compress = T, safe = T)

  # --

  outLib = list(PrecursorData = do.call('rbind', lapply(libList, add.idx, offsetIdx)),
                MsmsData = do.call('rbind',libList)[,grep(msmsCols,colnames(libList[[1]]))])
  print(paste("Saving processed library to: ",output,sep = ""))
  save(outLib, file = output, version = NULL, ascii = FALSE, compress = T, safe = T)
  gc(full = T)
  print("Processing completed!")

}
