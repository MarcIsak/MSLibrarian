#' Subsets a spectral library based on a selection of protein IDs
#' @param inputLib path to an input spectral library
#' @param outputLib desired path to a subset output spectral library
#' @param type means of selecting protein IDs. Possible values are 1) "calibrationlibrary" 2) "diann" or 3) custom character vector of protein IDs.
#' @param calibrationLib path to a calibration library object (*.RData). This argument is only used if type == "calibrationlibrary"
#' @param diaFiles character vector to DIA MS files in .raw format.
#' @param fdr numeric decimal value giving the protein group fdr cutoff to use for filtering (default = 0.2, 0 < fdr <= 1).
#' @param diannPath path to DIA-NN executable. Only used if type == "diann". If no path is given, the function will try to auto-detect the DIA-NN installation.
#' @param replace should the input library be overwritten/replaced (logical)? Default is FALSE
#' @export filter.proteins

filter.proteins <- function(inputLib = NULL, outputLib = NULL, type = NULL, calibrationLib = NULL, fdr = 0.2, diaFiles, diannPath = NULL, replace = F) {


  diann = F
  remove = F
  if(is.null(inputLib) | !file.exists(inputLib) | !str_detect(inputLib, ".tsv$|.csv$")) {
    stop("Argument 'inputLib' missing with the correct extensions (.tsv or .csv)")
  }
  if(is.null(type)) {
    stop("Arg - 'type' is missing.")
  } else if(type == "calibrationlibrary") {
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
        print(str_c("Unique protein groups in Calibration library: ", length(unique(calibLib@RawLib$ProteinId))))
        proteinIds = unique(unlist(str_split(calibLib@RawLib$ProteinId, ";")))
        proteinIds = str_remove(str_remove(proteinIds, "sp\\||tr\\|"), "\\|.*")
        str = "calib_lib"
      }
    }
  } else if(type == "diann") {
    diann = T
    str = "diann"
  } else if(class(type) == "character" & length(type) >= 1 & all(unique(nchar(type)) == 6 | unique(nchar(type)) == 10)) {
    print("Use a customized set of protein IDs to subset library...") # THIS IS NOT OPTIMAL. MUST BE MORE WAYS TO CHECK THAT THE INPUT IS UNIPROT IDENTIFIERS.
    proteinIds = type
    str = "custom"
  } else {
    stop("Arg - type. Invalid value provided. Run '?filter.proteins' in console to see valid values.")
  }
  if(is.null(outputLib)) {
    outputLib = file.path(dirname(inputLib), str_c(str_remove(basename(inputLib), ".tsv$|.csv$"),  "_protein_filter_", str,".",tools::file_ext(inputLib)))
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
    }
  }
  if(diann) {
    if(is.null(diannPath)) {
      print("Searching for DIA-NN installation...this may take a few seconds...")
      diannPath = system2("where", args = c("/r", "C:\\", "diann.exe"), stdout = T)
      if(length(diannPath) == 1) {
        print("Found DIA-NN executable...")
      } else {
        stop("Could not auto-detect the DIA-NN executable. Please add path to executable (diannPath argument).")
      }
    } else if(!(file.exists(diannPath) & str_detect(diannPath, "diann.exe"))) {
      stop("Specified DIA-NN executable does not exist...")
    } else {
      print(str_c("DIA-NN path: ", diannPath))
    }
    if(!dir.exists(file.path(dirname(inputLib), "diann"))) {
      dir.create(file.path(dirname(inputLib), "diann"))
      run.diann(diannPath = diannPath,
                diaFile = diaFiles,
                libFile = inputLib,
                output = file.path(dirname(inputLib), "diann"))
    }
    print("Loading DIA-NN report...")
    diannReport = read_tsv(file.path(dirname(inputLib), "diann", "report.tsv"),
                           col_types = cols_only(PG.Q.Value = col_double(),
                                                 Protein.Group = col_character()))
    if(max(diannReport$PG.Q.Value) <= fdr) {
      fdr = max(diannReport$PG.Q.Value)
      print(str_c("Maximum protein group FDR according to DIA-NN: ", fdr))
    }
    print(str_c("Extracting protein group ID at a FDR = ", fdr))
    proteinIds = diannReport$Protein.Group[diannReport$PG.Q.Value <= fdr]
    print(str_c("Unique protein groups in DIA-NN report after filtering: ", length(unique(proteinIds))))
    proteinIds = unique(unlist(str_split(proteinIds, ";")))
    proteinIds = str_remove(str_remove(proteinIds, "sp\\||tr\\|"), "\\|.*")
    outputLib = str_replace(outputLib, "_diann", str_c("_diann_", round(fdr, 2)))
  }
  protein.filter(inputLib = inputLib, outputLib = outputLib, proteinIds = proteinIds)
  if(remove) {
    file.remove(inputLib)
  }
  outputLib
}
