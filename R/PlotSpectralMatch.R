#' Visualizes spectral matches (experimental vs. predicted) in a mirror plot
#' @param msLib input MSLibrarian object
#' @param libType character defining the library type c("Consensus","Unfiltered"). Only spectral matches for the selected library type will be plotted.
#' @param specIdx numeric vector specifying the spectra indices to plot (if randSpec is FALSE).
#' @param trim logical indicating whether spectra should be mz-trimmed or not.
#' @param mzRange numeric vector specifying the range of mz-values to keep if trim = TRUE
#' @param tolType character specfifying if the tolerance should be "absolute" or "relative".
#' @param tol numeric value specifying the tolerance to use
#' @param precursorRm logical indicating if the precursor should be removed or not.
#' @param randIntensity logical determining if the Prosit relative intensities should be randomly permuted or not
#' @export plot.spectral.match

plot.spectral.match <- function(msLib, libType, specIdx, trim, mzRange, tolType, tol, precursorRm, randIntensity) {

  # #' @param randSpec logical indicating whether a random sample of spectral matches (exp vs. pred) should be plotted.
  # #' @param sampleSize integer defining the number of random spectral matches to plot (if randSpec is TRUE).


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
    relative= F
  } else if(tolType == "relative") {
    relative = T
  }
  if(trim == T) {
    trimRange = mzRange
    print("apply trim!")
  } else {
    trimRange = c(0,10^6) # It would be better to get the max mz value out from spectra objects
    print("no trim!")
  }

  plot.spectra <- function(idx,libType) {

    if(grepl("ox",slot(msLib@SpectrastLib,libType)@PrecursorData$FilterLib$ModifiedPeptideSequence[idx])) {
      modifications = c(C=57.021464,M=15.994915)
    } else {
      modifications = c(C=57.021464)
    }

    expSpec = trimMz(slot(msLib@SpectrastLib,libType)@Spectra[[slot(msLib@SpectrastLib,libType)@PrecursorData$FilterLib$msmsIdx[idx]]],
                     trimRange)
    predSpec = trimMz(msLib@PredLib@PrositLib$Spectra[[slot(msLib@SpectrastLib,libType)@PrecursorData$FilterLib$predIdx[idx]]],
                      trimRange)
    plot(expSpec, predSpec, relative = relative, tolerance = tol,
         sequences = c(slot(msLib@SpectrastLib,libType)@PrecursorData$FilterLib$PeptideSequence[idx],
                       msLib@PredLib@PrositLib$PrecursorData$StrippedPeptide[slot(msLib@SpectrastLib,libType)@PrecursorData$FilterLib$predIdx[idx]]),
         modifications = modifications)
    mtext(paste("Dotproduct",": ",round(compareSpectra(expSpec, predSpec,"dotproduct"),2)),side = 3, adj = 0.5)
    mtext(paste("Pearsons corr",": ",round(compareSpectra(expSpec, predSpec,"cor"),2)),side = 3, adj = 1)
    print(paste("Plotting index: ",as.character(idx),sep=""))

  }
  sapply(specIdx,plot.spectra,libType)
  # expSpec = trimMz(slot(msLib@SpectrastLib,libType)@Spectra[[slot(msLib@SpectrastLib,libType)@PrecursorData$FilterLib$msmsIdx[j]]],
  #                  trimRange)
  # predSpec = trimMz(msLib@PredLib@PrositLib$Spectra[[slot(msLib@SpectrastLib,libType)@PrecursorData$FilterLib$predIdx[j]]],
  #                   trimRange)

}
