#' Checks if DIA files are valid
#' @param diaFiles path to input DIA files

check.ms.files <- function(diaFiles) {
  print("Checking DIA files...")
  if(is.null(diaFiles)) {
    stop("Arg - diaFiles. No RAW files have been added")
  } else if(all(file.exists(diaFiles))) {
    diaFiles = diaFiles[grepl(".raw$|.d$", diaFiles)]
    if(length(diaFiles) == 0) {
      stop("No RAW files are supplied")
    } else if (any(basename(diaFiles) == diaFiles)) {
      stop("Arg - diaFiles. Full paths to RAW files must be supplied")
    }
  } else if(length(diaFiles) == 1){

    if(dir.exists(diaFiles) & length(list.files(diaFiles, pattern = ".raw$|.d$") != 0)) {
      print(str_c("Found raw files in folder:", diaFiles))
      diaFiles = list.files(diaFiles, pattern = ".raw$|.d$", full.names = T)
    } else {
      stop("Arg - diaFiles. Folder and/or raw files in folder do not exist")
    }
  } else {
    stop("Arg - diaFiles. Raw files or folder with raw files do not exist")
  }
  diaFiles
}
