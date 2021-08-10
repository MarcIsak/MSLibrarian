#' Converts MS RAW files to desired output format using MsConvert
#' @param msConvertPath absolute path to the MsConvert executable
#' @param rawFiles character vector with absolute paths to input raw files
#' @param format character giving the desired converted format
#' @param filter character specifing the filter to use during the conversion
#' @param output character giving the desired output folder for the converted files
#' @param threads integer giving the number of parallel conversion to run
#' @export run.msconvert

run.msconvert <- function(msConvertPath = NULL, rawFiles, format, filter, output, threads) {

  msconvert.cmd <- function(rawFile, msConvertPath, format, filter, output, verbose) {

    system2(msConvertPath, args = c(rawFile, format, filter, output, verbose),
            wait = T,
            stdout = file.path(str_remove(output,"-o "),
                               str_c(basename(rawFile), ".log"))) # Executes msconvert.exe with the options specified above.

  }
  if(is.null(msConvertPath)) {
    warning(str_c("Argument 'msConvertPath' is missing. Trying to auto-detect the MSconvert executable...may take a few seconds...", "\n"))
    msConvertPath = system2("where", args = c("/r", "C:\\", "msconvert.exe"), stdout = T)
    msConvertPath = msConvertPath[grep("ProteoWizard", msConvertPath)]
    if(length(msConvertPath) == 1 & identical(basename(msConvertPath), "msconvert.exe")) {
      print("Found MsConvert executable...")
    } else {
      stop("Could not auto-detect the MsConvert executable. Please add path to executable (msConvert argument).")
    }
  } else if(!(file.exists(msConvertPath) & str_detect(msConvertPath, "msconvert.exe"))) {
    stop("Specified MsConvert executable does not exist...")
  }
  if(threads > length(rawFiles)) {
    threads = length(rawFiles)
  }
  if(threads >= 8) {
    threads = 8
  }
  print(str_c("Setting threads to: ", threads))
  verbose = "-v"
  if (filter != "") {
    filter = str_c("--filter", shQuote(filter), sep = " ")
  }
  output = str_c("-o", output, sep = " ")
  cl = makeCluster(threads)
  clusterEvalQ(cl, {
    .libPaths(.libPaths())
    library(MSLibrarian)
  })
  parallel::parLapply(cl,
                      rawFiles,
                      msconvert.cmd,
                      msConvertPath,
                      format,
                      filter,
                      output,
                      verbose)
  stopCluster(cl)
  msConvertPath
}
