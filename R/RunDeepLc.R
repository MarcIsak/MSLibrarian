#' Runs DeepLC to predict retention times for a set of peptides
#' @param inputFolder path to folder with input files
#' @param splitCal split size for predictions
#' @param minDivider minimum dictionary divider
#' @param maxDivider maximum dictionary divider
#' @param deeplc path to deeplc_gui folder
#' @param threads number of threads to use for predictions
#' @export run.deeplc

run.deeplc <- function(inputFolder, splitCal = 50, minDivider = 10, maxDivider = 5000, deeplc, threads) {

  execute.deeplc <- function(predFile, calFile, threads, splitCal, minDivider, maxDivider) {
    system2(command = activate, args = c("deeplc_gui",
                                         "&",
                                         "deeplc",
                                         "--file_pred",
                                         predFile,
                                         "--file_cal",
                                         calFile,
                                         "--n_threads",
                                         threads,
                                         "--split_cal",
                                         splitCal,
                                         "--dict_divider",
                                         minDivider,
                                         "--dict_divider",
                                         maxDivider))
  }
  if(threads > 8) {
    threads = 8
  }
  if(is.null(deeplc)) {
    print("Searching for DeepLC installation...this may take a few seconds...")
    deeplcPath = system2(command = "where", args =  c("/r", "C:\\", "deeplc_gui.jar"),  stdout = T)
    activate  = file.path(dirname(deeplcPath), "Miniconda3/Scripts/activate.bat")
    if(length(deeplcPath) == 1 & file.exists(activate)) {
      print("Found DeepLC installation...")
    } else {
      stop("Could not auto-detect the DeepLC installation. Please add path to the DeepLC GUI folder (deeplc argument).")
    }
  } else if(!file.exists(file.path(deeplc, "Miniconda3/Scripts/activate.bat"))) {
    stop("Invalid path to DeepLC GUI folder...")
  } else {
    activate  = file.path(deeplcPath, "Miniconda3/Scripts/activate.bat")
  }

  lapply(c(file.path(inputFolder, "deeplc_lib.csv"),
           file.path(inputFolder, "deeplc_bench.csv")),
         execute.deeplc,
         file.path(inputFolder, "deeplc_calib.csv"),
         threads,
         splitCal,
         minDivider,
         maxDivider)

}
