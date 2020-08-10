#' Imports a processed MS2PIP library object into the current MSLibrarian object
#' @param msLib input MSLibrarian object
#' @param ms2pipLib file name of the library to be imported (incl. absolute path)
#' @export load.ms2pip.lib

load.ms2pip.lib <- function(msLib, ms2pipLib) {

  load(ms2pipLib)
  msLib@PredLib@Ms2pipLib = outLib
  msLib

}
