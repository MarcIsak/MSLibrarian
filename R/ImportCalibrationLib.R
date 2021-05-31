#' Creates a Spectrast Library object, and imports a Spectrast Library which is added to slot OriginalLib
#' @param library absolute path to the unfiltered library file (*.tsv)
#' @export import.calibration.lib

import.calibration.lib <- function(library) {

  calibLib = new("CalibrationLibrary", RawLib = as.data.frame(read_tsv(library, col_types = cols())))
  calibLib

}
