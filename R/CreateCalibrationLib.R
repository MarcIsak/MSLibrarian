#' Creates a calibration library from a set of DIA files and a FASTA file
#' @param diaFiles *Input DIA MS files in either RAW or mzXML format
#' @param fasta *FASTA file to use for database searching in Comet.
#' @param projectFolder *Name of project folder (absolute path). It must be unique.
#' @param msConvert path to msConvert executable. If no argument is passed, it will try to autodetect the executable on your system.
#' @param tppDir path to Transproteomic pipeline installation folder (for example "C:/TPP"). if no argument is passed, it will try to autodetect the installation folder.
#' @param openMsDir path to OpenMS installation folder.(such as "C:/Program Files/OpenMS-2.5.0/"). Will try to autodetect the installation folder if no argument is passed.
#' @param filter character vector with msCOnvert filters. Peakpickning on the MS1 level is set as default.
#' @param staggeredWindows logical stating if staggered windows were acquired for the input DIA files. (Default is FALSE).
#' @param diaUmpireParams *path to a DIA Umpire Signal extraction parameters file. These parameters are necessary for the creation of Pseudo-DDA MS/MS spectra.
#' @param cometParams *path to a Comet parameters file necessary for the database searching of pseudo-MS/MS spectra.
#' @param spectrastParams *path to a spectrast parameters file with parameters for creation of a spectral library.
#' @param libType Type of library to build. Possible values are 1) 'consensus' (default) and 2) 'best_replicate'
#' @param updateLib logical indicating if the created library should be updated. Will in this case re-run Spectrast and OpenSwathAssayGenerator to make a new library.
#' @param irt Can either be set to "biognosys_irt" or NULL (default). If set to "biognosys_irt", the calibration library retention times will be converted to iRT values.
#' @param threads Number of threads or parallel processes to use. Will by default use all available logical processors.
#' @export create.calibration.lib

create.calibration.lib <- function(diaFiles, fasta, projectFolder, msConvert = NULL, tppDir = NULL, openMsDir = NULL, filter = "peakPicking true 1-", staggeredWindows = F, diaUmpireParams = NULL, cometParams = NULL, spectrastParams = NULL, libType = "consensus", updateLib = FALSE, irt = NULL, threads = detectCores()) {


  find.params = function(params, paramsName, str) {
    if(is.null(params)) {
      warning(str_c("Argument '", paramsName, "' is missing. Searching for default parameter file...may take a few seconds...", "\n"))
      params = system2("where", args = c("/r", "C:\\", str), stdout = T)
      if(str_detect(params, str)) {
        print(str_c("Found default parameter file..."))
      } else {
        stop(str_c("Could not auto-detect the default parameter file. Please add path to this file using argument '", paramsName, "'"))
      }
    } else if(file.exists(params)){
      print(str_c("Using parameters file: ", params))
    } else {
      stop(str_c("Arg - '", paramsName,"'. File does not exist!"))
    }
    params
  }
  if(!updateLib) {
    if(!dir.exists(projectFolder)) {
      pass = dir.create(projectFolder)
      if(pass) {
        print(str_c("Creating project folder: ", projectFolder))
      } else {
        stop("Cannot create project folder")
      }
    } else {
      stop("Project folder already exists")
    }
  } else if (updateLib & dir.exists(projectFolder)) {
    print(str_c("Will update calibration library in project folder: ", projectFolder))
    prophetFolder = file.path(projectFolder, "prophet")
    outputFolder = file.path(projectFolder, "library")
  } else {
    print("Project folder does not exist. Cannot update the calibration library!")
  }
  if(is.null(tppDir)) {
    tppDir = dirname(system2("where", args = "PeptideProphetParser.exe",stdout = T))
    if(!dir.exists(tppDir[grep("TPP", tppDir)])) {
      stop("Cannot auto-dectect an installation of Trans-proteomic pipeline in your system. Please install or enter path to TPP installation folder!")
    } else {
      interact = system2("where", args = "InteractParser.exe",stdout = T)
      interact = interact[grep("TPP", interact)]
      comet = system2("where", args = "comet.exe",stdout = T)
      comet = comet[grep("TPP", comet)]
      peptideProphet = system2("where", args = "PeptideProphetParser.exe",stdout = T)
      peptideProphet = peptideProphet[grep("TPP", peptideProphet)]
      interProphet = system2("where", args = "InterProphetParser.exe",stdout = T)
      interProphet = interProphet[grep("TPP", interProphet)]
      spectrast = system2("where", args = "spectrast.exe",stdout = T)
      spectrast = spectrast[grep("TPP", spectrast)]
      if(!(length(interact) == 1 & length(peptideProphet) == 1 & length(interProphet) == 1 & length(comet) == 1 & length(spectrast) == 1)) {
        stop("Cannot find all TPP executables in TPP installation folder: InteractParser, PeptideProphet and InterProphet!")
      }
    }
  } else {
    comet = file.path(tppDir, "bin", "comet.exe")
    interact = file.path(tppDir, "bin", "InteractParser.exe")
    peptideProphet = file.path(tppDir, "bin", "PeptideProphetParser.exe")
    interProphet = file.path(tppDir, "bin", "InterProphetParser.exe")
    spectrast = file.path(tppDir, "bin", "spectrast.exe")
    if(!(file.exists(comet) & file.exists(interact) & file.exists(peptideProphet) & file.exists(interProphet) & file.exists(spectrast))) {
      stop("Cannot find all TPP executables in TPP installation folder: InteractParser, PeptideProphet and InterProphet!")
    }
  }
  if(is.null(openMsDir)) {
    openMsDir = dirname(system2("where", args = "OpenSwathAssayGenerator.exe",stdout = T))
    if(!dir.exists(openMsDir[grep("OpenMS", openMsDir)])) {
      stop("Cannot auto-dectect an installation of OpenMS in your system. Please install or enter path to TPP installation folder!")
    } else {
      oswGen = system2("where", args = "OpenSwathAssayGenerator.exe",stdout = T)
      oswGen = oswGen[grep("OpenMS", oswGen)]
      if(length(oswGen) != 1) {
        stop("Cannot find OpenSwathAssayGenerator.exe in OpenMS installation folder!")
      }
    }
  } else {
    oswGen = file.path(openMsDir, "bin","OpenSwathAssayGenerator.exe")

    if(!(file.exists(oswgen))) {
      stop("Cannot find OpenSwathAssayGenerator.exe in OpenMS installation folder!")
    }
  }
  if(!(libType == "consensus" | libType == "best_replicate")) {
    stop("Invalid arg - libType. Possible values are 1) 'consensus' or 2) 'best_replicate'")
  }
  cometParams = find.params(cometParams, paramsName = "cometParams", str = "comet_mslibrarian_default.params")
  diaUmpireParams = find.params(diaUmpireParams, paramsName = "diaUmpireParams", str = "diaumpire_se_mslibrarian_default.params")
  spectrastParams = find.params(spectrastParams, paramsName = "spectrastParams", str = "spectrast_create_mslibrarian_default.params")
  if(!updateLib) {
    if(all(str_detect(diaFiles, pattern = ".raw$"))) {

      if(staggeredWindows) {
        print("Staggered windows acquired for DIA runs...")
        filter = c(filter, "demultiplex optimization=overlap_only massError=10ppm")
      }
      tic()
      print("Converting the following DIA RAW files to mzXML format:")
      print(diaFiles)
      dir.create(file.path(projectFolder, "dia"))
      outputFolder = file.path(projectFolder, "dia")
      print(str_c("Creating output folder: ", outputFolder))
      msConvert = run.msconvert(msConvertPath = msConvert,
                                rawFiles = diaFiles, # msConvertPath = msConvert
                                format = "--mzXML",
                                filter = filter,
                                output = outputFolder,
                                threads = threads)
      diaFiles = list.files(outputFolder, pattern = ".mzXML", full.names = T)
      toc()
    }
    if(all(str_detect(diaFiles, pattern = ".mzXML$"))) {
      tic()
      print("Converting the following DIA files (mzXML) to pseudo-DDA MS/MS files (MGF): ")
      print(diaFiles)
      outputFolder = dir.create(file.path(projectFolder, "pseudo_dda"))
      outputFolder = file.path(projectFolder, "pseudo_dda")
      print(str_c("Creating output folder: ", outputFolder))
      make.pseudo.spectra(msConvertPath = msConvert,
                          mzmlFiles = diaFiles, # msConvertPath = msConvert
                          format = "--mgf",
                          diaUmpireParams = diaUmpireParams,
                          output = outputFolder,
                          threads = threads)
      toc()
    } else {
      stop("Incorrect MS files supplied. Make sure that all input DIA files are either in RAW or mzXML format.")
    }
    print("Preparing for Comet search...")
    tic()
    run.msconvert(msConvertPath = msConvert,
                  rawFiles = list.files(outputFolder, pattern = ".mgf$", full.names = T),
                  format = "--mzXML",
                  filter = "",
                  output = outputFolder,
                  threads = threads)
    toc()
    tic()
    print("Running Comet...")
    run.comet(cometPath = comet,
              msFiles = list.files(outputFolder, pattern = ".mzXML$", full.names = T), # file.path(tppDir, "bin", "comet.exe")
              fasta = fasta,
              cometParams = cometParams,
              threads = threads)
    toc()
    tic()
    print("Processing Comet results...")
    dir.create(file.path(projectFolder, "prophet"))
    prophetFolder = file.path(projectFolder, "prophet")
    print(str_c("Creating output folder: ", prophetFolder))
    run.prophets(interactPath = interact,
                 peptideProphetPath = peptideProphet,
                 interProphetPath = interProphet,
                 cometFiles = list.files(outputFolder, pattern = ".pep.xml$", full.names = T), # file.path(tppDir, "bin")
                 output = prophetFolder,
                 threads = threads)
    toc()
    tic()
    outputFolder = dir.create(file.path(projectFolder, "library"))
    outputFolder = file.path(projectFolder, "library")
    print(str_c("Creating output folder: ", outputFolder))

  }
  make.calibration.lib(spectrastPath = spectrast,
                       oswGenPath = oswGen,
                       prophetFile = file.path(prophetFolder, "interact_ipro.pep.xml"),
                       libType = libType,
                       spectrastParams = spectrastParams,
                       output = file.path(outputFolder, "calibration_lib.tsv"), # output = "D:/MS_files/MSLibrarian_DIA_test/brian_searle/calibration_lib/calibration_lib_dia.tsv",
                       irt = irt) # No irt file used...
  toc()
}
