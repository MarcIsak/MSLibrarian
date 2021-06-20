#' Runs DeepLC to predict retention times for a set of peptides
#' @param inputFolder path to folder with input files
#' @param splitCal split size for predictions
#' @param minDivider minimum dictionary divider
#' @param maxDivider maximum dictionary divider
#' @param deeplc path to deeplc_gui folder
#' @param threads number of threads to use for predictions
#' @export run.deeplc

run.deeplc <- function(inputFolder, splitCal = 50, minDivider = 10, maxDivider = 5000, deeplc, threads) {

  execute.deeplc <- function(predFile, activate, calFile, threads, splitCal, minDivider, maxDivider, gui) {

    if(gui) {
      pythonMod = c("deeplc_gui",
                    "&",
                    "python",
                    "-m",
                    "deeplc")
    } else {
      pythonMod = NULL
    }
    system2(command = activate, args = c(pythonMod,
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
                                         maxDivider),
            wait = T)
  }
  if(threads > 8) {
    threads = 8
  }
  if(is.null(deeplc)) {
    print("Searching for DeepLC GUI installation...this may take a few seconds...")
    deeplcPath = system2(command = "where", args =  c("/r", "C:\\", "deeplc_gui.jar"),  stdout = T)
    activate  = file.path(dirname(deeplcPath), "Miniconda3/Scripts/activate.bat")
    if(length(deeplcPath) == 1 & file.exists(activate)) {
      print("Found DeepLC GUI installation...")
      lapply(c(file.path(inputFolder, "deeplc_lib.csv"),
               file.path(inputFolder, "deeplc_bench.csv")),
             execute.deeplc,
             activate = activate,
             file.path(inputFolder, "deeplc_calib.csv"),
             threads,
             splitCal,
             minDivider,
             maxDivider,
             gui = T)
    } else {
      print("DeepLC GUI installation not found...")
      print("Searching for DeepLC executable instead...")
      deeplcPath = system2(command = "where", args =  c("/r", "C:\\", "deeplc.exe"),  stdout = T)
      if(length(deeplcPath) == 1 & file.exists(deeplcPath) & grepl("deeplc.exe$", deeplcPath)) {
        print("Found DeepLC executable...")
        activate = deeplcPath
        lapply(c(file.path(inputFolder, "deeplc_lib.csv"),
                 file.path(inputFolder, "deeplc_bench.csv")),
               execute.deeplc,
               activate = activate,
               file.path(inputFolder, "deeplc_calib.csv"),
               threads,
               splitCal,
               minDivider,
               maxDivider,
               gui = F)
      } else {
        stop("Could not auto-detect any DeepLC installation. Please add path to the DeepLC GUI folder (deeplc argument).")
      }
    }
  } else if(file.exists(file.path(deeplc, "Miniconda3/Scripts/activate.bat")) & grepl("deeplc_gui", deeplc)) {
    print("Found DeepLC GUI path...")
    lapply(c(file.path(inputFolder, "deeplc_lib.csv"),
             file.path(inputFolder, "deeplc_bench.csv")),
           execute.deeplc,
           activate = file.path(deeplcPath, "Miniconda3/Scripts/activate.bat"),
           file.path(inputFolder, "deeplc_calib.csv"),
           threads,
           splitCal,
           minDivider,
           maxDivider,
           gui = T)
    # stop("Invalid path to DeepLC GUI folder or DeepLC executable...")
  } else if (file.exists(deeplc) & grepl("deeplc.exe$", deeplc)) {
    print("Found path to DeepLC executable...")
    lapply(c(file.path(inputFolder, "deeplc_lib.csv"),
             file.path(inputFolder, "deeplc_bench.csv")),
           execute.deeplc,
           activate = deeplc,
           file.path(inputFolder, "deeplc_calib.csv"),
           threads,
           splitCal,
           minDivider,
           maxDivider,
           gui = F)

  } else {
    # activate  = file.path(deeplcPath, "Miniconda3/Scripts/activate.bat")
    stop("No valid path to DeepLC GUI folder or DeepLC executable...")
  }
}
