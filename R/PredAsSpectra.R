#' Creates a Spectra object (with Spectrum2 objects) for the selected predicted library
#' @param msLib input MSLibrarian object
#' @param selPred character specifying the predicted library to use c("Prosit", "MS2PIP")
#' @param threads number of processors to use for the parallel computation
#' @param outFile absolute path to output text file with console info from parallel workers
#' @export pred.as.spectra

pred.as.spectra <- function(msLib, selPred, threads, outFile) {

  if(selPred == "Prosit") {
    slotStr = "PrositLib"
  } else if (selPred == "MS2PIP") {
    slotStr = "Ms2pipLib"
  } else {
    print("This predicted library does not exist!")
  }

  idxRange = as.matrix(slot(msLib@PredLib, slotStr)$PrecursorData[,c("startIdx","endIdx")])

  cl <- makeCluster(threads, outfile = outFile)
  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths("C:/Users/marc/Dropbox (Human Neural Develop)/Marc/MISTR/R/win-library/3.5"))

  outSpectra = foreach(i = 1:nrow(slot(msLib@PredLib, slotStr)$PrecursorData), .packages = c("tcltk","MSnbase")) %dopar% {

    if(!exists("pb")) {
      pb = tkProgressBar("Parallel task",min = 1, max = nrow(slot(msLib@PredLib, slotStr)$PrecursorData))
    }

    setTkProgressBar(pb, i)
    new("Spectrum2",
        precursorMz = slot(msLib@PredLib, slotStr)$PrecursorData$PrecursorMz[i],
        mz = slot(msLib@PredLib, slotStr)$MsmsData$FragmentMz[min(idxRange[i,]):max(idxRange[i,])],
        rt = slot(msLib@PredLib, slotStr)$PrecursorData$iRT[i],
        intensity = slot(msLib@PredLib, slotStr)$MsmsData$RelativeIntensity[min(idxRange[i,]):max(idxRange[i,])],
        centroided = T,
        smoothed = F,
        msLevel = as.integer(2))
  }
  stopCluster(cl)
  slot(msLib@PredLib, slotStr)$Spectra = Spectra(outSpectra)
  msLib
    # Can also add the iRT information directly....
}


