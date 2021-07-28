#' Runs the MSFragger database search tool on pseudo-DDA MS/MS spectra
#' @param msfragger path to MSFragger executable
#' @param fragParams path to MSFragger parameters file
#' @param msFiles input pseudo-dda MS/MS files
#' @param fasta input fasta file
#' @param decoyDb path to DecoyDatabase executable (part of OpenMS)

run.msfragger <- function(msfragger, fragParams, msFiles, fasta, decoyDb) {

  create.decoy.db <- function(decoyDb, fasta, decoyFasta) {

    system2(decoyDb, args = c("-in",
                                 fasta,
                                 "-out",
                                 decoyFasta,
                                 "-decoy_string",
                                 "DECOY_",
                                 "-decoy_string_position",
                                 "prefix",
                                 "-type",
                                 "protein",
                                 "-method",
                                 "reverse",
                                 "-enzyme",
                                 "Trypsin"))

  }
  run.msfragger <- function(fragCmd, fragParams, msFiles) {

    system2(command = "java", args = c("-Xmx8g",
                                       "-jar",
                                       fragCmd,
                                       fragParams,
                                       msFiles))

  }
  decoyFasta = str_replace(fasta, "\\.fasta$", "_decoy_concatenated.fasta")
  print("Adding decoys to FASTA...")
  create.decoy.db(decoyDb = decoyDb,
                  fasta = fasta,
                  decoyFasta = decoyFasta)
  print("Preparing MSFragger parameters file")
  params = as.data.frame(str_replace(t(read_delim(file = fragParams,
                                                  delim = "\t")),
                                     pattern = "database_name = test.fasta",
                                     replacement = str_c("database_name = ",
                                                         decoyFasta)))
  fragParams = str_replace(fragParams,
                           pattern = "\\.params",
                           replacement = "_tmp.params")
  colnames(params) = "# MSFragger"
  write_delim(x = params, file = fragParams, delim = "\t", quote_escape = "none")
  run.msfragger(fragCmd = msfragger,
                fragParams = fragParams,
                msFiles = msFiles)


}
