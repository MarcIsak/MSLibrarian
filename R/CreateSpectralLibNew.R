#' Creates a spectral library with intensities and retention times from Prosit
#' @param calibrationLib path to processed calibration library (RData format)
#' @param projectFolder path to project folder containing a processed calibration library.
#' @param fasta path to protein sequence FASTA file from Uniprot
#' @param msFile path to representative pseudo ms/ms file used to create the calibration library (could be extracted from the calibration library file later...)
#' @param outputLib path to desired output library
#' @param format output library format c("spectronaut", "openswath")
#' @param ceMode Collision energy optimization setting. Possible values are c("length_only", "length_charge", "charge", scalar(fixed collision energy))
#' @param saveObj Determines whether a MSLibrarian object will be saved. Possible values are "TRUE" or "FALSE" (default)
#' @export create.spectral.lib.new

create.spectral.lib.new <- function(calibrationLib = NULL, projectFolder = NULL, msFile = NULL, fasta, outputLib = NULL, format = "spectronaut", ceMode = "length_charge", threads = detectCores(), saveObj = F) {

  if(!is.null(projectFolder) & is.null(calibrationLib)) {
    calibrationLib = file.path(projectFolder, "library", "calibration_lib.RData")
    if(!dir.exists(projectFolder)) {
      stop("Project folder does not exist!")
    } else if(!file.exists(file.path(projectFolder, "library", "calibration_lib.RData"))) {
      stop(str_c("Could not find the library file: ", calibrationLib))
    } else {
      print(str_c("Found the library file: ", calibrationLib, " in project folder: ",projectFolder))
    }
    if(is.null(msFile)) {
      msFile = list.files(file.path(projectFolder, "dia"), pattern = ".mzXML$", full.names = T)[1]
      if(!file.exists(msFile)) {
        stop(str_c("Could not find any mzXML file(s) in: ", dirname(msFile)))
      } else {
        print(str_c("Found mzXML file(s) in ", dirname(msFile)))
      }
    } else if (!file.exists(msFile) | !str_detect(msFile, ".mzXML$")) {
      stop("Argument msFile: MS file does not exist with the correct extension (.mzXML)")
    } else {
      print(str_c("Found MS file: ", msFile))
    }
  } else if (is.null(calibrationLib) | !file.exists(calibrationLib) | !str_detect(calibrationLib, ".RData$")) {
    stop("Could not find a processed calibration library with the correct extension (.RData)")
  } else {
    print(str_c("Found processed calibration library file: ", calibrationLib))
  }
  if(is.null(msFile) | !file.exists(msFile) | !str_detect(msFile, ".mzXML")) {
    stop(str_c("Argument msFile - MS file does not exist with the correct extension (.mzXML)"))
  } else {
    print(str_c("Using MS file: ", msFile))
  }
  if(is.null(outputLib)) {
    outputLib = file.path(dirname(calibrationLib),
                          str_c(format(Sys.time(),
                                       "%b_%d_%H_%M_%S_%Y"),
                                "_mslibrarian_ce_",
                                ceMode, "_irt_prosit"))
    print(str_c("No output library file specified! Generating default: ", outputLib))
  } else if (!dir.exists(dirname(outputLib))) {
    stop("Specified output library folder does not exists!")
  }
  if(threads > detectCores() | threads > 8) {
    threads = 8
  }
  tic()
  msLib = read.fasta(file = fasta) # The FASTA FILE CAN ONLY HANDLE UNIPROT/SWISS-PROT DATABASES AT THE MOMENT...
  toc()
  tic()
  msLib = digest.proteins(msLib,  # See if the mass calculation step can be speed up. It is the time-limiting step + very slow...
                          rowStr = "Accession",
                          enzyme = "trypsin",
                          maxMissed = 0,
                          carbamidomethyl = T,
                          threads = threads)
  toc()
  tic()
  msexp = readMSData(msFile, mode = "onDisk")
  msLib = filter.peptides.tmp(msLib,
                              AALengthRange = c(7, 30),
                              rejectAAs = c("U", "O", "B", "J", "X", "Z"), # Can be taken from the calibration library.
                              mzRange = range(fData(msexp)[fData(msexp)$msLevel == 1, c("lowMZ", "highMZ")]),
                              chargeRange = c(2,3)) # A parameter should be set for this later, but it will cause a bug at the moment.
  toc()
  tic()
  load(calibrationLib)
  msLib = map.db.entries(msLib,
                         predDb = calibLib@MetaData$PredictionDb,
                         chargeRange = c(2,3)) # A parameter should be set for this later, but it will cause a bug at the moment.
  toc()
  tic()
  msLib = get.precursors.new(msLib,
                             mzRange = range(fData(msexp)[fData(msexp)$msLevel == 1, c("lowMZ", "highMZ")]),
                             chargeRange = c(2,3), # A parameter should be set for this later, but it will cause a bug at the moment.
                             matchDb = T,
                             threads = threads)
  toc()
  tic()
  msLib = make.spectra.lib(msLib,
                           bestCe = get.best.ce(ceMat = calibLib@Comparisons),
                           predDb = calibLib@MetaData$PredictionDb,
                           ceMode = ceMode)
  toc()
  tic()
  outputLib = export.spectra.lib(msLib, # It must return the output lib argument...
                                 calibrationLib = calibrationLib,
                                 format = format,
                                 outputLib = outputLib,
                                 batchSize = 10000000)
  if(saveObj) {
    print(str_c("Saving MSLibrarian object to: ", str_replace(outputLib, pattern = ".tsv$|.csv$", replacement = ".RData")))
    save(msLib, file = str_replace(outputLib, pattern = ".tsv$|.csv$", replacement = ".RData"),
         version = NULL,
         ascii = FALSE,
         compress = T,
         safe = T)
  }
  rm(msLib)
  gc()
  toc()

}
