#' Adds MS/MS information from Spectral library file into a MSnbase Spectra object (slot: Spectra)
#' @param calibrationLib input SpectrastLibrary object
#' @export make.msms2

make.msms2 <- function(calibrationLib) {

  add.msms.data <- function(prec, spectrum, msmsData) {

    idx = as.integer(prec["msmsIdx"])
    spectrum[[idx]]@precursorMz = as.numeric(prec["PrecursorMz"])
    spectrum[[idx]]@precursorCharge = as.integer(prec["PrecursorCharge"])
    spectrum[[idx]]@rt = as.numeric(prec["NormalizedRetentionTime"])
    spectrum[[idx]]@mz = msmsData[as.numeric(prec["startIdx"]):as.numeric(prec["endIdx"]),"ProductMz"]
    spectrum[[idx]]@intensity = msmsData[as.numeric(prec["startIdx"]):as.numeric(prec["endIdx"]),"LibraryIntensity"]
    spectrum[[idx]]@intensity = spectrum[[idx]]@intensity[order(spectrum[[idx]]@mz)]
    spectrum[[idx]]@mz = spectrum[[idx]]@mz[order(spectrum[[idx]]@mz)]
    spectrum[[idx]]@peaksCount = length(spectrum[[idx]]@intensity)
    spectrum[[idx]]@centroided = T
    spectrum[[idx]]

  }
  precData = calibrationLib@PrecursorData$Original
  msmsData = calibrationLib@RawLib[,c("ProductMz","LibraryIntensity")]
  print("Adding precursor and MS/MS information to Spectrum2 objects...")
  spectra = replicate(new("Spectrum2"), n = nrow(precData))
  calibrationLib@Spectra = apply(precData, 1, add.msms.data, spectra, msmsData)
  calibrationLib
}
