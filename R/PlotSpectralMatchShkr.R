#' Visualizes spectral matches (experimental vs. predicted) in a mirror plot
#' @param pepShkr input PeptideShaker object
#' @param prositLib input Prosit spectral library file
#' @param indx numeric index speciying the PSM to plot to corresponding predicted spectra.
#' @param trim logical indicating if spectra should be mz-trimmed or not.
#' @param mzRange numeric vector specifying the range of mz-values to include in the matching.
#' @param precursorRm logical indicating if the precursor should be removed or not.
#' @param mode character specifying whether "dotProduct" or "cor" should be annotated in the plot
#' @param toleranceType character specfifying if the tolerance should be "absolute" or "relative".
#' @param tolerance numeric value specifying the tolerance to use
#' @param randIntensity logical determining if the Prosit relative intensities should be randomly permuted or not
#' @export plot.spectral.match.shkr

plot.spectral.match.shkr <- function(pepShkr,prositLib,indx,trim,mzRange,precursorRm,mode,toleranceType,tolerance,randIntensity) {

  rm.precursor <- function(x) {

    k = T
    while(k == T) {

      f = mz(x)[which(intensity(x) == max(intensity(x)))]

      if(precursorMz(x) > (f-1) & precursorMz(x) < (f+1)) {
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

  if(trimMz == F) {
    mzRange = c(0,10^6)
  }

  predSpec = trimMz(prositLib$Spectrum2[[indx]], mzRange)
  expSpec = trimMz(pepShkr@ms2Spectra[[pepShkr@psmsPass$specLib$`spectrum title indx`[indx]]],mzRange)

  if(precursorRm == T) {
    expSpec = rm.precursor(expSpec)
  }
  if(randIntensity == T){
    predSpec = do.random.intensity(predSpec)
  }
  if(toleranceType == "absolute") {
    relative = F
  } else {
    relative = T
  }

  plot(expSpec,predSpec,relative = relative,tolerance = tolerance, sequences= c(pepShkr@psmsPass$specLib$peptide_ref[indx],
                                     prositLib$Precursor$StrippedPeptide[pepShkr@psmsPass$specLib$prositIndx[indx]]))
  mtext(paste(mode,": ",round(compareSpectra(expSpec, predSpec,mode),2)),side = 3, adj = 1)

}
