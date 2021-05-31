#' Creates input csv files for retention time predictions in DeepLC
#' @param inputLib path to input spectral library in CSV-format (Spectronaut) which should have its retention time changed.
#' @param calibLib path to calibration library as created by the MSLibrarian package (.RData)
#' @param outputFolder path to output folder where DeepLC input files will be added.
#' @param nCal numeric value (nCal > 0 & nCal <= 1), giving the proportion of peptides with known retention times to use for calibration.
#' @export make.deeplc.inputs

make.deeplc.inputs <- function(inputLib, calibLib, outputFolder, nCal) {

  add.mod <- function(seq) {

    if(str_detect(seq, "C")) {
      paste(as.character(unique(unlist(str_locate_all(seq,"C")))),
            "Carbamidomethyl",
            sep = "|",
            collapse = "|")
    } else {
      ""
    }
  }
  resultFolder = file.path(outputFolder, str_replace(basename(inputLib), ".tsv|.csv", str_c("_deeplc_", nCal)))
  print(str_c("Creating results folder: ", resultFolder))
  dir.create(resultFolder)

  if(str_detect(inputLib, ".csv$")) {
    print(str_c("Importing spectronaut library : ", inputLib))
    specLib = read_csv(inputLib, col_types = cols_only(iRT = col_double(),
                                                       ModifiedPeptide = col_character(),
                                                       PrecursorCharge = col_integer(),
                                                       StrippedPeptide = col_character()))
    delim = ","

  } else if(str_detect(inputLib, ".tsv$")) {
    print(str_c("Importing OpenSwath library: ", inputLib))
    specLib = read_tsv(inputLib, col_types = cols_only(NormalizedRetentionTime = col_double(),
                                                       ModifiedPeptideSequence = col_character(),
                                                       PrecursorCharge = col_integer(),
                                                       PeptideSequence = col_character()))
    delim = "\t"
    colnames(specLib)[grep("NormalizedRetentionTime", colnames(specLib))] = "iRT"
    colnames(specLib)[grep("ModifiedPeptideSequence", colnames(specLib))] = "ModifiedPeptide"
    colnames(specLib)[grep("PeptideSequence", colnames(specLib))] = "StrippedPeptide"
  }
  precursors = str_c(specLib$ModifiedPeptide, "_", specLib$PrecursorCharge)
  libSeqs = specLib$StrippedPeptide[which(!duplicated(precursors))] # Would be better in the future to use ModifiedPeptideColumn instead.
  print(str_c("Found: ", length(libSeqs), " precursors in the input library..."))
  print("Importing the calibration library...")
  rm(specLib)
  gc()
  benchData = calibLib@PrecursorData$FilterLib[,c("PeptideSequence", "NormalizedRetentionTime")]
  nSeq = ceiling(nrow(benchData)*nCal)
  print(str_c("Sampling ", nSeq, " peptide sequences ( ",round(nCal*100,2)," % ) ",
              "with known retention times to use for calibration of DeepLC predictions..."))
  calibIdx = sample(seq(1,nrow(benchData)), nSeq)

  print("Writes csv file for prediction of retention times for library peptides")
  libPeptides = write_csv(data.frame(seq = libSeqs,
                                     modifications = sapply(libSeqs, add.mod)),
                          file = str_c(resultFolder, "/", "deeplc_lib.csv"))
  print("Writes csv file for prediction of retention times for benchmark peptides")
  benchPeptides = write_csv(data.frame(seq = benchData$PeptideSequence,
                                       modifications = sapply(benchData$PeptideSequence, add.mod)),
                            file = str_c(resultFolder, "/","deeplc_bench.csv"))
  print("Writes csv file for calibration of the retention time predictions")
  calibPeptides = write_csv(data.frame(seq = benchData$PeptideSequence[calibIdx],
                                       modifications = sapply(benchData$PeptideSequence[calibIdx], add.mod),
                                       tr = benchData$NormalizedRetentionTime[calibIdx]),
                            file = str_c(resultFolder, "/","deeplc_calib.csv"))
  resultFolder

}
