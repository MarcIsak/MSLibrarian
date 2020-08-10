#' Loads a Spectrast Library object into the MSLibrarian object
#' @param msLib input MSLibrarian object
#' @param file absolute path to Spectrast Library object (saved as RData)
#' @export load.spectrast.lib

load.spectrast.lib <- function(msLib, file) {

  load(file)
  msLib@SpectrastLib = spectrastLib
  msLib

}
