#' Loads a processed Prosit Library file into a MSLibrarian object
#' @param msLib input MSLibrarian object
#' @param file absolute path to processed Prosit library file
#' @export load.prosit.lib

load.prosit.lib <- function(msLib, file) {

  load(file)
  msLib@PredLib@PrositLib = outLib
  msLib

}
