#' Make part of prosit library file into Spectrum2 objects contained in a Spectra object
#' @param msLib input MSLibrarian object
#' @param threads number of processors to use for the parallel computation
#' @param outFile absolute path to output text file with console info from parallel workers
#' @export prosit.as.spectra

prosit.as.spectra <- function(msLib, threads, outFile) {

  idxRange = as.matrix(msLib@PredLib@PrositLib$PrecursorData[,c("startIdx","endIdx")])

  cl <- makeCluster(threads, outfile = outFile)
  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths("C:/Users/marc/Dropbox (Human Neural Develop)/Marc/MISTR/R/win-library/3.5"))

  outSpectra = foreach(i = 1:nrow(msLib@PredLib@PrositLib$PrecursorData), .packages = c("tcltk","MSnbase")) %dopar% {

    if(!exists("pb")) {
      pb = tkProgressBar("Parallel task",min = 1, max = nrow(msLib@PredLib@PrositLib$PrecursorData))
    }

    setTkProgressBar(pb, i)
    new("Spectrum2",
        precursorMz = msLib@PredLib@PrositLib$PrecursorData$PrecursorMz[i],
        mz = msLib@PredLib@PrositLib$MsmsData$FragmentMz[min(idxRange[i,]):max(idxRange[i,])],
        rt = msLib@PredLib@PrositLib$PrecursorData$iRT[i],
        intensity = msLib@PredLib@PrositLib$MsmsData$RelativeIntensity[min(idxRange[i,]):max(idxRange[i,])],
        centroided = T,
        smoothed = F,
        msLevel = as.integer(2))
  }
  stopCluster(cl)
  msLib@PredLib@PrositLib$Spectra = Spectra(outSpectra)
  msLib
    # Can also add the iRT information directly....
}


