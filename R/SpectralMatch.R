#' Compares experimental spectra (Spectrast) with corresponding predicted spectra (for example Prosit)
#' @param msLib input MSLibrarian object
#' @param tolType A character specifying if matching peak tolerance should be "relative" or "absolute"
#' @param tol numeric specifying the matching tolerance
#' @param trim logical indicating if trim should be performed or not
#' @param mzRange numeric vector indicating mz values to be removed
#' @param precursorRm logical indicating whether to remove intense precursors to improve spectral matching
#' @param randIntensity logical indicating if random permutation of intensities in Prosit Spectrum2 obj. should be performed.
#' @param listName character specifying the list name for accessing the metrics of the spectra comparison
#' @param threads number of processors to use for the parallel processing
#' @param outFile absolute path to the output text file with console info from the parallel workers
#' @export spectral.match

spectral.match <- function(msLib,tolType,tol,trim,mzRange,precursorRm,randIntensity, listName, threads, outFile) {

  rm.precursor <- function(x) {

    k = T
    while(k == T) {

      maxInt = mz(x)[which(intensity(x) == max(intensity(x)))]

      if(precursorMz(x) > (maxInt-1) & precursorMz(x) < (maxInt+1)) {
        q = intensity(x)
        q[which(intensity(x) == max(intensity(x)))] = 0
        x@intensity  = q
        print("Intense Precursor!")
      } else {
        k = F
      }
    }
    x
  }
  do.random.intensity <- function(x) {
    x@intensity = sample(intensity(x),replace = F)
    x
  }

  if(tolType =="absolute") {
    tt= F
  } else if(tolType == "relative") {
    tt = T
  }
  if(trim == T) {
    trimRange = mzRange
    print("apply trim!")
  } else {
    trimRange = c(0,10^6) # It would be better to get the max mz value out from spectra objects
    print("no trim!")
  }

  # slotStr = c("Consensus", "Unfiltered")
  slotStr = "Unfiltered"

  for (i in 1:length(slotStr)) {

    cl <- makeCluster(threads, outfile = outFile)
    registerDoParallel(cl)
    clusterEvalQ(cl, .libPaths("C:/Users/marc/Dropbox (Human Neural Develop)/Marc/MISTR/R/win-library/3.5"))
    # nrow(slot(msLib@SpectrastLib,slotStr[i])@PrecursorData$FilterLib)
    matches = foreach(j = 1:nrow(slot(msLib@SpectrastLib,slotStr[i])@PrecursorData$FilterLib), .packages = c("tcltk","MSnbase"), .combine = 'rbind') %dopar% {

      if(!exists("pb")) {
        pb = tkProgressBar("Parallel task",min = 1, max = nrow(slot(msLib@SpectrastLib,slotStr[i])@PrecursorData$FilterLib))
      }

      setTkProgressBar(pb, j)

      if(!is.na(slot(msLib@SpectrastLib,slotStr[i])@PrecursorData$FilterLib$predIdx[j])) {

        expSpec = trimMz(slot(msLib@SpectrastLib,slotStr[i])@Spectra[[slot(msLib@SpectrastLib,slotStr[i])@PrecursorData$FilterLib$msmsIdx[j]]],
                         trimRange)
        predSpec = trimMz(msLib@PredLib@PrositLib$Spectra[[slot(msLib@SpectrastLib,slotStr[i])@PrecursorData$FilterLib$predIdx[j]]],
                          trimRange)
        if(precursorRm == T){
          expSpec = rm.precursor(expSpec)
        }
        if(randIntensity == T) {
          predSpec = do.random.intensity(predSpec)
        }
        c(compareSpectra(expSpec, predSpec, relative = tt, tolerance = tol, "common"),
          compareSpectra(expSpec, predSpec,"dotproduct"),
          compareSpectra(expSpec, predSpec,"cor"))
      } else {

        c(NA,NA,NA)
      }
    }
    stopCluster(cl)
    colnames(matches) = c("commonPeaks","dotProduct","pearsonsCor")
    gc()
    msLib@Comparisons[[slotStr[i]]][[listName]] = matches
  }
  gc()
  msLib
}
