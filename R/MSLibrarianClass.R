#' Class definition of a MSLibrarian object
#'

setClass("CalibrationLibrary", slots = c(RawLib = "data.frame",
                                         PrecursorData = "list",
                                         Spectra = "list",
                                         Comparisons = "data.frame",
                                         MetaData = "list"))

setClass("Sequence", slots = c(Proteins = "data.frame",
                               Peptides = "list",
                               Precursors = "list",
                               LibraryMatch = "data.frame"))

setClass("PredictedLibrary", slots = c(PrositLib = "list",
                                       Ms2pipLib = "list",
                                       DeepMassPrism = "list",
                                       DeepDIA = "list",
                                       pDeep2 = "list",
                                       DeepRTPlus = "list"))
setClass("MSLibrarian", slots = c(Sequences = "Sequence",
                                  PredLib = "PredictedLibrary",
                                  CalibLib= "CalibrationLibrary"))
# setClass("MSLibrarian", slots = c(Sequences = "Sequence",
#                             PredLib = "PredictedLibrary",
#                             SpectrastLib = "SpectrastLibrary",
#                             Comparisons = "list"))




# The SpectrastLibrary class is set in package SpectrastLib, added as a dependency to this package.






