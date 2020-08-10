#' Creates a Spectra object (with Spectrum2 objects) for the selected predicted library
#' @param msLib input MSLibrarian object
#' @param selPred character specifying the predicted library to use c("Prosit", "MS2PIP")
#' @param threads number of processors to use for the parallel computation
#' @param outFile absolute path to output text file with console info from parallel workers
#' @export pred.as.spectra.tmp

pred.as.spectra.tmp <- function(msLib, selPred, threads, outFile) {

  if(selPred == "Prosit") {
    slotStr = "PrositLib"
    selPred = "prositIdx"
  } else if (selPred == "MS2PIP") {
    slotStr = "Ms2pipLib"
    selPred = "ms2pipIdx"
  } else {
    print("This predicted library does not exist!")
  }

  idxRange = as.matrix(slot(msLib@PredLib, slotStr)$PrecursorData[,c("startIdx","endIdx")])

  cl <- makeCluster(threads, outfile = outFile)
  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths("C:/Users/marc/Dropbox (Human Neural Develop)/Marc/MISTR/R/win-library/3.5"))
  predIdx = msLib@SpectrastLib@Unfiltered@PrecursorData$FilterLib[,selPred]

  outSpectra = foreach(i = 1:length(predIdx), .packages = c("tcltk","MSnbase")) %dopar% {

    if(!exists("pb")) {
      pb = tkProgressBar("Parallel task",min = 1, max = length(predIdx))
    }

    setTkProgressBar(pb, i)
    if(!is.na(predIdx[i])) {
      tmp = new("Spectrum2",
            precursorMz = slot(msLib@PredLib, slotStr)$PrecursorData$PrecursorMz[predIdx[i]],
            mz = slot(msLib@PredLib, slotStr)$MsmsData$FragmentMz[min(idxRange[predIdx[i],]):max(idxRange[predIdx[i],])],
            rt = slot(msLib@PredLib, slotStr)$PrecursorData$iRT[predIdx[i]],
            intensity = slot(msLib@PredLib, slotStr)$MsmsData$RelativeIntensity[min(idxRange[predIdx[i],]):max(idxRange[predIdx[i],])],
            centroided = T,
            smoothed = F,
            msLevel = as.integer(2))
    } else {
      new("Spectrum2")
    }
  }
  stopCluster(cl)
  slot(msLib@PredLib, slotStr)$Spectra = Spectra(outSpectra)
  #outSpectra
  msLib
    # Can also add the iRT information directly....
}


