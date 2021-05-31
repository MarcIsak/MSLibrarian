#' Replaces the retention times for an input spectral library
#' @param inputLib input spectral library that will have its retention times replaced.
#' @param outputLib path to output library
#' @param calibrationLib calibration library object used to extract calibration peptides for DeepLC
#' @param outputFolder desired path to output folder for DeepLC results
#' @param nCal decimal number giving the proportion of calibration peptides in calibrationLib to use (default = 0.33)
#' @param deeplc path to the DeepLC GUI installation folder (will try to auto-detect it by default)
#' @param replace should the input library be replaced(logical)? Default is FALSE
#' @export replace.rts

replace.rts <- function(inputLib = NULL, outputLib = NULL, calibrationLib = NULL, outputFolder = NULL, nCal = 0.33, deeplc = NULL, replace = F, threads = detectCores()) {

  if(is.null(inputLib) | !file.exists(inputLib) | !str_detect(inputLib, ".tsv$|.csv$")) {
    stop("Argument 'inputLib' missing with the correct extensions (.tsv or .csv)")
  }
  if(is.null(calibrationLib)) {
    stop("Arg - calibrationLib is missing.")
  } else if(!file.exists(calibrationLib) | !str_detect(calibrationLib, ".RData$")) {
    stop("File does not exist and/or does not have the correct extension (must be .RData)")
  } else {
    load(calibrationLib)
    if(class(calibLib) != "CalibrationLibrary") {
      stop("Arg - calibrationLib. This object is not of class CalibrationLibrary. Please import the correct RData file.")
    } else {
      print(str_c("Importing calibration library: ", calibrationLib))
    }
  }
  if(is.null(outputFolder)) {
    outputFolder = dirname(inputLib)
    print(str_c("Output folder set to: ", outputFolder))
  } else {
    if(dir.exists(outputFolder)) {
      print(str_c("Output folder: ", outputFolder))
    } else {
      stop("Output folder does not exist!")
    }
  }
  if(is.null(outputLib)) {
      if(str_detect(basename(inputLib), "irt_prosit")) {
        outputLib = file.path(dirname(inputLib), str_replace(basename(inputLib), "irt_prosit", str_c("rt_deeplc_", nCal)))
      } else {
        outputLib = file.path(dirname(inputLib), str_c(str_remove(basename(inputLib), ".tsv$|.csv"), "_rt_deeplc_", nCal, ".", tools::file_ext(inputLib)))
      }
  } else if(identical(inputLib, outputLib) & !replace) {
    stop("Overwriting the input library is forbidden when argument replace is set to false")
  }
  if(replace) {
    print("Arg - replace = TRUE")
    if(identical(inputLib, outputLib)) {
      print("Will overwrite the input library after replacing retention times.")
      outputLib = inputLib
    } else {
      print("Will remove input library after processing...")
      remove = T
    }
  }
  if(!identical(file_ext(inputLib), file_ext(outputLib))) {
    stop("Output library must have the same extension as the input library")
  } else {
    if(!dir.exists(dirname(outputLib))) {
      stop("Output directory does not exist...")
    } else {
      print(str_c("Writing filtered output library to: ", outputLib))
    }
  }
  deepLcFolder = make.deeplc.inputs(inputLib = inputLib, calibLib = calibLib, outputFolder = outputFolder, nCal = nCal)
  run.deeplc(inputFolder = deepLcFolder, deeplc = deeplc, threads = threads)
  add.deeplc.rt(inputLib = inputLib, outputLib = outputLib, calibLib = calibLib, nCal, deepLcFolder = deepLcFolder)
  if(remove) {
    file.remove(inputLib)
  }

}
