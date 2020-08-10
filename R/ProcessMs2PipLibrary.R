#' Processes a MS2PIP predicted library and saves it as a RData object
#' @param zipfilePred The absolute path to the zipped intensity predictions output file from MS2PIP (in spectronaut.csv format)
#' @param zipfileRT The absolute path to the zipped retention time predictions output file from DeepLC (in spectronaut.csv format)
#' @param output Output file name for the RData object (including absolute path)
#' @export process.ms2pip.lib

process.ms2pip.lib <- function(zipfilePred, zipfileRT, output) {

  get.csv <- function(zipf, zipPath) {

    counter <<- counter + 1
    print(paste("Loading library file: ", as.character(counter), sep = ""))
    y = read.csv(unz(zipPath, zipf), quote = "", stringsAsFactors = F, header = T)

  }
  ms2pip.list <- function(x) {
    tmpLib[tmpLib[,"PrecursorCharge"] == x,]
  }
  get.idx <- function(x) {
    cols = "StrippedPeptide|ModifiedPeptide|PrecursorMz|PrecursorCharge"
    #startIdx = which(!duplicated(x[,"ModifiedPeptide"])) # This does not work with methionine residues I think...
    startIdx = which(!duplicated(x[,"ProteinId"])) # This works!!!
    endIdx = c(startIdx[2:length(startIdx)]-1,nrow(x))
    cbind(startIdx,endIdx,x[startIdx,grep(cols,colnames(x))])
  }
  get.set.idx <- function(x,setIdx) {
    setRow = sum(setIdx[1:i])
    x$startIdx = x$startIdx + setRow
    x$endIdx = x$endIdx + setRow
    i <<- i + 1
    x
  }

  print("Initiates processing of the MS2PIP library")
  counter <<- 0
  fileNames = unzip(zipfile = zipfilePred, list = TRUE)$Name
  tmpLib = do.call('rbind', lapply(fileNames, get.csv, zipfilePred))
  tmpLib = tmpLib[order(proteinId = as.numeric(str_replace(tmpLib$ProteinId, "peptide",""))),]
  libList = lapply(unique(tmpLib[,"PrecursorCharge"]), ms2pip.list)
  rm(tmpLib)
  gc()
  print("Indexing of spectral library...")
  print("Removal of redundant entries...")
  setIdx = c(0,as.numeric(lapply(libList,nrow)))
  spectraCols = "RelativeIntensity|FragmentMz|FragmentType|FragmentNumber|FragmentCharge"
  i <<- 1
  print("Merging library files into a single library...")
  outLib = list(PrecursorData = do.call('rbind',lapply(lapply(libList,get.idx),get.set.idx, setIdx)),
                MsmsData = do.call('rbind',libList)[,grep(spectraCols,colnames(libList[[1]]))])
  print("Adding predicted iRT values to the library...")
  rtPreds = read.csv(unz(zipfileRT, unzip(zipfileRT, list = T)$Name),quote = "", stringsAsFactors = F, header = T)

  if(all(rtPreds$seq == outLib$PrecursorData$StrippedPeptide)) {
    print(all(rtPreds$seq == outLib$PrecursorData$StrippedPeptide)) # remove after
    outLib$PrecursorData$iRT = rtPreds$predicted_tr
  } else {
    #print(length(rtPreds$seq)) # Remove after
    #print(length(outLib$PrecursorData$StrippedPeptide)) # remove after
    print("Peptide sequences do not match between library and iRT predictions --> iRT predictions are not added...")
  }

  print(paste("Saving processed library to: ",output,sep = ""))
  save(outLib, file = output, version = NULL, ascii = FALSE, compress = T, safe = T)
  gc()
  print("Processing completed!")
  list(irt = rtPreds, lib = outLib$PrecursorData) # Remove after

}

