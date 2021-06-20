#' Create a calibration library in OpenSwath format (TSV)
#' @param spectrastPath absolute path to the Trans-proteomic pipeline binary folder
#' @param oswGenPath absolute path to the OpenMS binary folder
#' @param prophetFile absolute path to the iProphet results file
#' @param libType Type of spectral library to build. Possible values are 1) 'Consensus' and 2) 'Best replicate'"
#' @param irt input file with irt values associated to every peptide sequence
#' @param spectrastParams input Spectrast parameters file
#' @param output Desired output file name of OpenSwath library
#' @export make.calibration.lib

make.calibration.lib <- function(spectrastPath, oswGenPath, prophetFile, libType, irt, spectrastParams, output) {

  if(libType == "consensus") {
    libMode = "-cAC"
    print("Consensus spectra chosen for building spectral library!")
  } else if(libType == "best_replicate") {
    libMode = "-cAB"
    print("Best replicate spectra chosen for building spectral library!")
  } else {
    print("Invalid library type. Possible values are 1) 'Consensus' and 2) 'Best replicate'")
  }
  if(is.null(irt)) {
    print("No iRT file provided...")
    irt = NULL
  } else if(irt == "biognosys_irt"){
    print("Library retention times will be converted to Biognosys iRT values...")
    write_delim(data.frame(sequence = c("LGGNEQVTR", "GAGSSEPVTGLDAK", "VEATFGVDESNAK", "YILAGVENSK",
                                        "TPVISGGPYEYR", "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK",
                                        "GTFIIDPGGVIR", "GTFIIDPAAVIR", "LFLQFGAQGSPFLK"),
                           iRT = c(-24.916114, 0, 12.3893748888888, 19.7879106666667,
                                   28.7145812222221, 33.381243, 42.2638884444446, 54.621042, 70.5187413333333,
                                   87.2332223333333, 100.002821666667)),
                file = file.path(dirname(prophetFile), "irt_biognosys.txt"),
                delim = "\t",
                col_names = F)
    irt = str_c("-c_IRT", file.path(dirname(prophetFile), "irt_biognosys.txt"))
  } else {
    print("invalid argument irt...will ignore this parameter...")
    irt = NULL
  }
  #if(!file.exists(str_replace(prophetFile, ".pep.xml$", ".splib"))) {
  print("Creating Spectrast library in SPLIB format")
  system2(spectrastPath, # str_c(tppPath, "spectrast.exe")
          args = c(str_c("-cN", str_remove(prophetFile, pattern = ".pep.xml")),
                   irt,
                   str_c("-cF", spectrastParams),
                   prophetFile))
  #}
  splibFile = str_replace(prophetFile, ".pep.xml$", ".splib")
  print(str_c("Found Spectrast library: ", splibFile))
  print(str_c("Library build action: ", libType))
  system2(spectrastPath, # str_c(tppPath, "spectrast.exe")
          args = c("-cM",
                   libMode,
                   str_c("-cN", str_replace(prophetFile, ".pep.xml$", str_c("_", libType))),
                   irt,
                   str_c("-cF", spectrastParams),
                   splibFile))
  inLib = str_replace(prophetFile, ".pep.xml", str_c("_", libType, ".mrm"))
  system2(oswGenPath,
          args = c("-max_transitions 100",
                   "-product_lower_mz_limit 200.0",
                   str_c("-in ", inLib),
                   str_c("-out ", output)))

}








# assayGenPath = "C:/Program Files/OpenMS-2.5.0/bin/OpenSwathAssayGenerator.exe"
# mrmInput = paste("-in", "D:/MS_files/MSLibrarian_DIA_test/comet_files/interact_ipro.mrm")
# tsvOutput = paste("-out", "D:/MS_files/MSLibrarian_DIA_test/comet_files/interact_ipro.tsv")
# maxTransitions = "-max_transitions 100"
# productLowerMzLimit = "-product_lower_mz_limit 200.0"
#
# system2(assayGenPath, args = c(maxTransitions, productLowerMzLimit, mrmInput, tsvOutput))
