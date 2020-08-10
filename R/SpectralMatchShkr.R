#' Returns dot-product, pearsons correlation and common peaks for matches between PSMs in PepShkr object and predicted spectra from Prosit.
#' @param pepShkr input pepShkr object
#' @param prositLib input spectral library file from Prosit
#' @param toleranceType A character specifying if matching peak tolerance should be "relative" or "absolute"
#' @param tol numeric specifying the matching tolerance
#' @param trim logical indicating if trim should be performed or not
#' @param mzRange numeric vector indicating mz values to be removed
#' @param listName Name of list to add results to
#' @param precursorRm logical indicating whether to remove intense precursors to improve spectral matching
#' @param randIntensity logical indicating if random permutation of intensities in Prosit Spectrum2 obj. should be performed.
#' @export spectral.match.shkr

spectral.match.shkr <- function(pepShkr, prositLib, toleranceType, tol, trim, mzRange,listName,precursorRm,randIntensity) {

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

  matchCol = c("commonPeaks","dotProduct","pearsonsCor")
  matches = matrix(NA,nrow(pepShkr@psmsPass$specLib),length(matchCol))
  colnames(matches) = matchCol

  if(toleranceType =="absolute") {
    tt= F
  } else if(tolerance == "relative") {
    tt = T
  }
  if(trim == T) {
    trimRange = mzRange
    print("apply trim!")
  } else {
    trimRange = c(0,10^6) # It would be better to get the max mz value out from spectra objects
    print("no trim!")
  }

  for (i in 1:nrow(pepShkr@psmsPass$specLib)) {

    if(!is.na(pepShkr@psmsPass$specLib$prositIndx[i])) {

      tmp = trimMz(pepShkr@ms2Spectra[[pepShkr@psmsPass$specLib$`spectrum title indx`[i]]],trimRange)
      prositTmp = trimMz(prositLib$Spectrum2[[i]],trimRange)

      if(precursorRm == T){
        tmp = rm.precursor(tmp)
      }
      if(randIntensity == T) {
        prositTmp = do.random.intensity(prositTmp)
      }
      matches[i,"commonPeaks"] = compareSpectra(tmp, prositTmp,relative = tt, tolerance = tol, "common")
      matches[i,"dotProduct"] = compareSpectra(tmp, prositTmp,"dotproduct")
      matches[i,"pearsonsCor"] = compareSpectra(tmp, prositTmp,"cor")
    }
    print(i)
  }
  pepShkr@psmsPass[[listName]] = as.data.frame(matches)
  pepShkr

}
