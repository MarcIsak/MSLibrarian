#' Creates pseudo DDA spectra f
#' @param msConvertPath absolute path to the MsConvert executable
#' @param mzmlFiles character vector with absolute paths to input raw files
#' @param format character giving the desired converted format
#' @param diaUmpireParams character specifying the absolute path to the diaUmpire SE parameters file
#' @param output character giving the desired output folder for the converted files
#' @param threads integer giving the number of parallel conversion to run
#' @export make.pseudo.spectra

make.pseudo.spectra <- function(msConvertPath, mzmlFiles, format, diaUmpireParams, output, threads) {

  run.msconvert <- function(mzmlFile, msConvertPath, format, filter, output, verbose) {

    system2(msConvertPath, args = c(mzmlFile, format, filter, output, verbose),
            wait = T,
            stdout = file.path(str_remove(output,"-o "),
                               str_c(basename(mzmlFile), ".log")),
            stderr = file.path(str_remove(output,"-o "),
                               str_c(basename(mzmlFile), ".log"))) # Executes msconvert.exe with the options specified above.

  }
  if(is.null(msConvertPath)) {
    msConvertPath = system2("where", args = c("/r", shQuote("C:\\Program Files"), "msconvert.exe"), stdout = T)
    msConvertPath = msConvertPath[grep("ProteoWizard", msConvertPath)]
    if(length(msConvertPath) == 1) {
      print("Found MsConvert executable...")
    } else {
      stop("Could not auto-detect the MsConvert executable. Please add path to executable (msConvert argument).")
    }
  } else if(!(file.exists(msConvertPath) & str_detect(msConvertPath, "msconvert.exe"))) {
    stop("Specified MsConvert executable does not exist...")
  }
  filter = str_c("--filter ", shQuote(str_c("diaUmpire params=", diaUmpireParams))) # A filter for creating pseudo-ms/ms spectra
  output = str_c("-o ", output) # Folder path for converted files
  verbose = "-v" # Detailed description of the conversion process
  if (threads == 1) {
    print("Parallel conversion disabled...")
    system2(msConvertPath, args = c(mzmlFiles, format, filter, output, verbose))
  } else if (threads >= length(mzmlFiles)) {
    threads = length(mzmlFiles)
    if(threads >= 8) {
      threads = 8
    }
    print(str_c("Parallel conversion enabled. Number of parallel processes: ", threads, "..."))
    cl = makeCluster(threads)
    clusterEvalQ(cl, {
      .libPaths(.libPaths())
      library(MSLibrarian)
    })
    parallel::parLapply(cl, mzmlFiles, run.msconvert, msConvertPath, format, filter, output, verbose)
    stopCluster(cl)
  } else {
    print("Incorrect value for parameter called threads...terminating code execution")
    stop()
  }
  gc()
  #hell
}
