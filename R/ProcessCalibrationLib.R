#' Creates and processes a calibration library object
#' @param specLib path to spectral library (in TSV-format, as processed by make.calibration.lib). Not necessary if the project Folder argument is added
#' @param projectFolder Path to project folder. If not specified it is necessary to add the specLib argument.
#' @param predictionDb path to SQLite database with Prosit predictions.
#' @param rt character specifying the input library retention time format. Possible values: c("seconds", "iRT")
#' @param outLib desired path to the output calibration library (.RData format)
#' @param threads number of  threads to use
#' @export process.calibration.lib

process.calibration.lib <- function(specLib = NULL, projectFolder = NULL, predictionDb, rt = NULL, outLib = str_replace(specLib, ".tsv", ".RData"), threads = detectCores()) {

  if(!is.null(projectFolder)) {
    if(dir.exists(projectFolder) & file.exists(file.path(projectFolder, "library", "calibration_lib.tsv"))) {
      specLib = file.path(projectFolder, "library", "calibration_lib.tsv")
      outLib = str_replace(specLib, ".tsv$", ".RData")
      print(str_c("Processing library: ", specLib))
    } else {
      stop("Could not find the project folder with the spectral library file (.TSV) as created by create.calibration.lib()")
    }
  } else if (is.null(specLib) | !str_detect(specLib, ".tsv$") | !file.exists(specLib)) {
    stop("Add specLib argument to a calibration library with the correct extension (.TSV)")
  } else {
    print(str_c("Processing library: ", specLib))
  }
  if(!dir.exists(dirname(outLib))) {
    stop("Output directory does not exist")
  } else if(!str_detect(outLib, ".RData$")) {
    outLib = str_c(outLib, ".RData")
  }
  if(!file.exists(predictionDb) | !str_detect(predictionDb, ".sqlite$")){
    stop("Argument missing for SQLite database with prosit predictions.")
  }
  if(is.null(rt)) {
    stop("No rt argument added. Please add.")
  } else if(rt == "sec") {
    k = T
    print("Library has retention times in seconds")

    } else if(rt == "iRT") {
      print("Library has indexed retention times (iRT)")
      k = F
    } else {
      stop("Invalid rt argument. Possible values are: 'sec' or 'iRT'")
    }
  if(threads > detectCores() | threads > 8) {
    threads = 8
  }

  calibLib = import.calibration.lib(library = specLib)

  calibLib = get.lib.precursors2(calibrationLib = calibLib)

  calibLib = filter.lib.precursors2(calibrationLib = calibLib,
                                    exUniMod = c(1),
                                    precAARange = c(7,30),
                                    chargeRange = c(2,3),
                                    exAAs = c("U", "O", "B", "J", "X", "Z"))
  if(k) {
    calibLib@PrecursorData$FilterLib$NormalizedRetentionTime =
      calibLib@PrecursorData$FilterLib$NormalizedRetentionTime/60
  }

  calibLib = make.msms2(calibrationLib = calibLib) # Could use threads here to speed it up...

  calibLib@Comparisons = run.spectral.match2(spectra = calibLib@Spectra,
                                             precursors = calibLib@PrecursorData$FilterLib[,c("PeptideSequence", "PrecursorCharge", "msmsIdx")],
                                             predictionDb = predictionDb,
                                             threads = threads)
  calibLib@MetaData$PredictionDb = predictionDb
  print(str_c("Writing processed Calibration Library to: ", outLib))
  save(calibLib, file = outLib,
       version = NULL,
       ascii = FALSE,
       compress = T,
       safe = T)
  ce.bench.plot(ceMat = calibLib@Comparisons, outLib = outLib)
  gc()
  #calibLib = rep.unimod.acc(calibLib)

}
