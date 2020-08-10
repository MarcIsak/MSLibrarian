#' Loads a predicted spectral library and extracts the relevant information such as retention time, relative intensities of fragment ions etc.
#' @param zipfile absoute path to the zip file containing the spectral library file
#' @param filename name of the file to be extracted from the zip file
#' @param output absolute path to processed library file, as a RData object
#' @export process.prosit.lib

process.prosit.lib <- function(zipfile,filename, output) {

  prosit.list <- function(x) {
    myPrositLib[myPrositLib[,"PrecursorCharge"] == x,]
  }
  get.set.idx <- function(x,setIdx) {
    setRow = sum(setIdx[1:i])
    x$startIdx = x$startIdx + setRow
    x$endIdx = x$endIdx + setRow
    i <<- i + 1
    x
  }
  get.idx <- function(x) {
    cols = "StrippedPeptide|LabeledPeptide|PrecursorMz|iRT|PrecursorCharge"
    startIdx = which(!duplicated(x[,"LabeledPeptide"]))
    endIdx = c(startIdx[2:length(startIdx)]-1,nrow(x))
    cbind(startIdx,endIdx,x[startIdx,grep(cols,colnames(x))])
  }
  print("Importing Prosit spectral library...")
  myPrositLib = read.csv(unz(zipfile,filename), quote="", stringsAsFactors=FALSE)
  print("Sorting spectral library by precursor charges...")
  libList = lapply(unique(myPrositLib$PrecursorCharge),prosit.list)
  rm(myPrositLib)
  gc()
  print("Indexing of spectral library...")
  print("Removal of redundant entries...")
  i <<- 1
  setIdx = c(0,as.numeric(lapply(libList,nrow)))
  spectraCols = "RelativeIntensity|FragmentMz|FragmentType|FragmentNumber|FragmentCharge"
  outLib = list(PrecursorData = do.call('rbind',lapply(lapply(libList,get.idx),get.set.idx, setIdx)),
              MsmsData = do.call('rbind',libList)[,grep(spectraCols,colnames(libList[[1]]))])
  print(paste("Saving processed library to: ",output,sep = ""))
  save(outLib, file = output, version = NULL, ascii = FALSE, compress = T, safe = T)
  gc()
  print("Processing completed!")

}
