#' Modifies and filters a spectral library on protein, peptide, or transition level. Also allows for changing library retention times.
#' @param projectFolder path to a MSLibrarian project folder
#' @param inputLib path to an input spectral library. If there is only one library, it will be selected automatically if arg - projectFolder is supplied.
#' @param outputLib desired path to a modified spectral library. A default name is generated if not specified.
#' @param calibrationLib path to a Calibration Library (*.RData format)
#' @param mods character vector specifying the modifications to apply. c("protein", "peptide", "transitions", "rt")
#' @param diaFiles path to DIA MS raw files used to build the calibration library
#' @param protMod type of protein modification to apply. Possible values are 1) "calibrationlibrary", 2) "diann" or 3) a character vector of Uniprot identfiers
#' @param protFdr decimal number between 0 - 1 giving the protein group FDR to use for protein ID subsetting of a library.
#' @param diannPath path to DIA-NN executable. If not supplied, the function will try to auto-detect the executable on your system.
#' @param topPept integer giving the maximum number of high responding peptides (as given by PREGO), to allow for each protein group.
#' @param pregoPath path to PREGO executable. If not supplied, the function will try to auto-detect the executable on your system.
#' @param cutoffTrans decimal number between 0 - 1 giving the minimum relative intensity that a transition can have (0.05 is default).
#' @param topTrans maximum number of high intensity transitions to allow for each precursor in the library. (Is ignored/NULL by default)
#' @param minTrans minimum number of transitions to allow for each precursor in a library.
#' @param nCal decimal number giving the proportion of calibration library entries to use for RT calibration in DeepLC. (0.25 by default)
#' @param deepLcPath path to DeepLC GUI folder. It will try to auto-detect the folder if not supplied.
#' @export mod.spectral.lib

mod.spectral.lib <- function(projectFolder = NULL, inputLib = NULL, outputLib = NULL, calibrationLib = NULL, mods, diaFiles, protMod, protFdr = 0.05, diannPath = NULL, topPept = NULL, pregoPath = NULL, cutoffTrans = NULL, topTrans = NULL, minTrans = NULL, nCal = 0.25, deepLcPath = NULL) {

  if(!is.null(projectFolder) & is.null(calibrationLib)) {
    calibrationLib = file.path(projectFolder, "library", "calibration_lib.RData")
    if(!dir.exists(projectFolder)) {
      stop("Project folder does not exist!")
    } else if(!file.exists(calibrationLib)) {
      stop(str_c("Could not find the library file: ", calibrationLib))
    } else {
      print(str_c("Found the library file: ", calibrationLib, " in project folder: ",projectFolder))
    }
  } else if(is.null(calibrationLib)) {
    stop("Arg - calibrationLib is missing.")
  } else if(!file.exists(calibrationLib) | !str_detect(calibrationLib, ".RData$")) {
    stop("File does not exist and/or does not have the correct extension (must be .RData)")
  } else {
    load(calibrationLib)
    if(class(calibLib) != "CalibrationLibrary") {
      stop("Arg - calibrationLib. This object is not of class CalibrationLibrary. Please import the correct RData file.")
    } else {
      print(str_c("Found processed calibration library file: ", calibrationLib))
      rm(calibrationLib)
      gc()
    }
  }
  if(is.null(inputLib) & !is.null(projectFolder)) {
    inputLib = list.files(file.path(projectFolder, "library"), pattern = ".tsv$|.csv$", full.names = T)
    inputLib = inputLib[!grepl("calibration_lib", inputLib)]
    if(length(inputLib) > 1) {
      stop(str_c("More than 1 spectral library found in: ", file.path(projectFolder, "library"),". Please add only one of these libraries with arg - inputLib."))
    } else if(length(inputLib) == 0) {
      stop(str_c("No spectral library not found in: ", file.path(projectFolder, "library"), ". Use arg - 'inputLib' to specify a spectral library to modify."))
    } else {
      print(str_c("Input spectral library found: ", inputLib, " in project folder: ", file.path(projectFolder, "library")))
    }
  } else if(is.null(inputLib)) {
    stop("Arg - inputLib. Please enter a path to a spectral library to modify.")
  } else if(!file.exists(inputLib) | !str_detect(inputLib, ".tsv$|.csv$")) {
    stop("Arg - inputLib. File does not exist and/or does not have the correct extension (.tsv or .csv)")
  } else {
    print(str_c("Input spectral library: ", inputLib))
  }
  if(!all(grepl("^protein$|^peptide$|^transition$|^rt$", mods))) {
    stop(str_c("Arg - mods: ", mods[!grepl("protein|peptide|transition|rt", mods)], " is not a valid value..."))
  }
  replace = F
  if(any(str_detect(mods, "^protein$"))) {
    print("Applying protein ID subsetting...")
    inputLib = filter.proteins(inputLib = inputLib,
                               type = protMod,
                               diaFiles = diaFiles,
                               fdr = protFdr,
                               calibrationLib = calibrationLib,
                               diannPath = diannPath,
                               replace = replace)
    replace = T
  }
  if(any(str_detect(mods, "^peptide$"))) {
    print("Applying peptide subsetting...")
    inputLib = filter.prego(inputLib = inputLib,
                            top = topPept,
                            pregoPath = pregoPath,
                            replace = replace)
    replace = T
  }
  if(any(str_detect(mods, "^transition$"))) {
    print("Applying transition filtering...")
    inputLib = filter.transitions(inputLib = inputLib,
                                  topTrans = topTrans,
                                  cutoffTrans = cutoffTrans,
                                  minTrans = minTrans,
                                  replace = replace)
    replace = T
  }
  if(any(str_detect(mods, "^rt$"))) {
    print("Replacement of retention times enabled...")
    replace.rts(inputLib = inputLib,
                calibrationLib = calibrationLib,
                nCal = nCal,
                replace = replace,
                deeplc = deepLcPath)
    replace = T
  }

}

