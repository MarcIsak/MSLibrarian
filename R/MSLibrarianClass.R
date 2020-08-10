#' Class definition of a MSLibrarian object
#'

setClass("MSLibrarian", slots = c(Sequences = "Sequence",
                            PredLib = "PredictedLibrary",
                            SpectrastLib = "SpectrastLibrary",
                            Comparisons = "list"))

setClass("Sequence", slots = c(Proteins = "data.frame",
                               Digestion = "list",
                               Peptides = "list",
                               FilterPeptides = "list",
                               Precursors = "data.frame",
                               LibraryMatch = "data.frame"))

setClass("PredictedLibrary", slots = c(PrositLib = "list",
                                       Ms2pipLib = "list",
                                       DeepMassPrism = "list",
                                       DeepDIA = "list",
                                       pDeep2 = "list",
                                       DeepRTPlus = "list"))

# The SpectrastLibrary class is set in package SpectrastLib, added as a dependency to this package






