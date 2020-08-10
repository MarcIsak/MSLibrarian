#' Make part of prosit library file into Spectrum2 objects contained in a Spectra object
#' @param pepShkr input PepShkr object
#' @param prositLib input prosit library file
#' @export prosit.as.spectrum2

prosit.as.spectrum2 <- function(pepShkr, prositLib) {

  prositLib[["Spectrum2"]] = Spectra()

  for (i in 1:nrow(pepShkr@psmsPass$specLib)) {

    if(!is.na(pepShkr@psmsPass$specLib$prositIndx[i])) {

    indxRange = as.numeric(c(prositLib$Precursor[pepShkr@psmsPass$specLib$prositIndx[i],c("startIndx","endIndx")])) # Extract indices for prosit intensities and mz values
    prositLib$Spectrum2[[i]] = new("Spectrum2",precursorMz = prositLib$Precursor$PrecursorMz[pepShkr@psmsPass$specLib$prositIndx[i]],
              mz = prositLib$Spectra$FragmentMz[seq(min(indxRange),max(indxRange))],
              rt = prositLib$Precursor$iRT[pepShkr@psmsPass$specLib$prositIndx[i]],
              intensity = prositLib$Spectra$RelativeIntensity[seq(min(indxRange),max(indxRange))],
              centroided = T,smoothed = F,msLevel = as.integer(2))
    # Can also add the iRT information directly....
    } else {

      prositLib$Spectrum2[[i]] = new("Spectrum2") # Empty Spectrum2 object
    }
    print(i)
  }
  prositLib
}


