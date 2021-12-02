#' Runs a DIA-NN analysis of the selected DIA-file and input library in TSV format
#' @param diannPath Absolute path to the DIA-NN executable
#' @param diaFile Absolute path to the input DIA raw file
#' @param libFile Absolute path to the input spectral library file (as given by the make.osw.lib function)
#' @param output Desired file name for DIA-NN reports (incl. absolute path)
#' @param matrices A logical determining if quantitative matrices for proteins and precursors should be outputted (Default=TRUE).
#' @param matQ numeric in the range c(0,1) setting the q-value to use for filtering matrices.
#' @param threads Integer specifying the number of threads for the analysis
#' @param fasta Absolute path to a FASTA file for matching of peptides sequences to proteins
#' @param relaxProtInf Logical determining if relaxed protein inference should be enabled. Default = FALSE
#' @param enzyme Character specifying the enzyme to use for in-silico digestion of sequences in the FASTA file
#' @param separateMassAcc should mass accuracy be determined for each individual MS file (logical)
#' @param massAcc Numeric setting the desired MS1 accuracy
#' @param reportLibInfo Logical determining if extra info on fragment ions and precursors should be saved to main report (default TRUE)
#' @param robustLC Specfies whether Quantification should be run in Robust LC mode. Possible values are "high_accuracy" or "high_precision". Default = NULL.
#' @param carbamidomethyl_C Carbamidomethyl on Cysteine residues as fixed modification. Default = T.
#' @param nativeLib Logical setting if the spectral library is a native library or not.
#' @param matchBetweenRuns logical determining if match-between-runs (MBR) should be enabled. (Default = FALSE)
#' @export run.diann

run.diann <- function(diannPath, diaFile, libFile, output = NULL, matrices = T, matQ = 0.01, threads = detectCores(), fasta = NULL, enzyme, separateMassAcc = F, massAcc = NULL, reportLibInfo = T, robustLC = NULL, relaxProtInf = F, carbamidomethyl_C = T, nativeLib = F, matchBetweenRuns = F) {

  if(separateMassAcc) {
    separateMassAcc = "--individual-mass-acc"
  } else {
    separateMassAcc = NULL
  }
  if(!is.null(massAcc)) {
    if(is.numeric(massAcc)) {
      print(str_c("Sets MS1 accuracy to: ", massAcc))
      massAcc = str_c("--mass-acc-ms1 ", massAcc)
    } else {
      stop("Arg - massAcc. Invalid value!")
    }
  }
  if(is.null(robustLC)) {
    peakCenter = NULL
  } else if(robustLC == "high_accuracy") {
    peakCenter = "--peak-center"
    print("Using Robust LC (high accuracy mode)")
  } else if(robustLC == "high_precision"){
    peakCenter = str_c("--peak-center", " --no-ifs-removal")
    print("Using Robust LC (high precision mode)")
  } else {
    stop("Arg - robustLC. Invalid string provided. Valid strings are 'high_accuracy' or 'high_precision'")
  }
  if(reportLibInfo) {
    print("Add extra info on fragment ions and precursors to main report...")
    reportLibInfo = "--report-lib-info"
  } else {
    reportLibInfo = NULL
  }
  if(carbamidomethyl_C) {
    carbamidomethyl = "--unimod4"
  } else {
    carbamidomethyl = NULL
  }
  if(matrices) {
    print("Will write quantitative protein and precursor matrices...")
    matrices = "--matrices"
    if(is.numeric(matQ) & matQ > 0 & matQ < 1) {
      print(str_c("Matrices will be filtered at Q = ", matQ))
      matQ = str_c("--matrix-qvalue ", matQ)
    } else {
      stop("Invalid q-value set for matrix filtering!")
    }
  } else {
    matrices = NULL
  }
  if(!is.null(fasta)) {
    fasta = str_c("--fasta ", fasta)
    noProtInf = NULL
    if(enzyme == "trypsin") {
      enzyme = "--cut K*,R*,!*P"
    } else {
      stop("Incorrect enzyme specified...code execution terminated")
    }
  } else {
    print("No FASTA...")
    noProtInf = "--no-prot-inf"
    enzyme = NULL
  }
  if(relaxProtInf) {
    print("The same protein cannot be in two separate groups...")
    relaxProtInf = "--relaxed-prot-inf"
  } else {
    relaxProtInf = NULL
  }
  # if(matchBetweenRuns) {
  #   print("Match-between-runs enabled!")
  #   mbr = "--mbr"
  # } else {
  #   mbr = NULL
  # }
  if(str_detect(libFile, ".tsv$") & !nativeLib) {
    # Should have some function that checks that the input library has the correct data...
    if(sum(str_detect(colnames(read_tsv(libFile, n_max = 1, col_types = cols())), "ModifiedPeptideSequence|LibraryIntensity|NormalizedRetentionTime")) == 3) {
      print("OpenSwath format selected")
      headers = paste(c("ModifiedPeptideSequence", "PrecursorCharge", "PrecursorMz", "NormalizedRetentionTime",
                        "ProductMz", "LibraryIntensity", "UniprotId", "*", "*", "*", "Decoy", "ProductCharge",
                        "FragmentType", "FragmentSeriesNumber"), collapse = ",")
      headers = str_c("--library-headers ", headers)
    } else {
      stop("TSV-file is not an OpenSwath library!")
    }

  } else if(str_detect(libFile, ".tsv$") & nativeLib) {
    colsOne = c("PrecursorMz","ProductMz", "Annotation","ProteinId", "GeneName", "PeptideSequence", "ModifiedPeptideSequence","PrecursorCharge",
             "LibraryIntensity", "NormalizedRetentionTime", "PrecursorIonMobility", "FragmentType", "FragmentCharge", "FragmentSeriesNumber",
             "FragmentLossType")
    colsTwo = c("FileName", "PrecursorMz", "ProductMz", "Tr_recalibrated", "IonMobility", "transition_name", "LibraryIntensity",
                "transition_group_id", "decoy", "PeptideSequence", "Proteotypic", "QValue", "PGQValue", "Ms1ProfileCorr", "ProteinGroup",
                "ProteinName", "Genes", "FullUniModPeptideName", "ModifiedPeptide", "PrecursorCharge", "PeptideGroupLabel", "UniprotID",
                "NTerm", "CTerm", "FragmentType", "FragmentCharge", "FragmentSeriesNumber", "FragmentLossType", "ExcludeFromAssay")
    if(identical(colnames(read_tsv(libFile, n_max = 1, col_types = cols())), colsOne) | identical(colnames(read_tsv(libFile, n_max = 1, col_types = cols())), colsTwo)) {
      print("Native TSV-library")
      headers = NULL
    } else {
      stop("Native library has incorrect column names...")
    }
  } else if(str_detect(libFile, ".csv$")) {
    # Should have some function that checks that the input library has the correct data...
    if(sum(str_detect(colnames(read_csv(libFile, n_max = 1, col_types = cols())), "ModifiedPeptide|RelativeIntensity|iRT")) == 3) {

      headers = paste(c("ModifiedPeptide", "PrecursorCharge", "PrecursorMz", "iRT",
                        "FragmentMz", "RelativeIntensity", "UniprotId", "*", "*", "*", "*", "FragmentCharge",
                        "FragmentType", "FragmentNumber"), collapse = ",")
      headers = str_c("--library-headers ", headers)
    } else {
      stop("CSV-file is not a Spectronaut library!")
    }
    print("Spectronaut format selected")
  } else if(str_detect(libFile, ".speclib")) {
    print(str_c("DIA-NN internal spectral library format selected: ", libFile))
    headers = NULL
  } else {
    stop("Arg - libFile. No library provided!")
  }
  print("Start analysis in DIA-NN!")
  system2(diannPath, args = c(paste("--f", diaFile),
                              paste("--lib", libFile),
                              paste("--threads", as.character(threads)),
                              fasta,
                              enzyme,
                              headers,
                              paste("--out", file.path(output, "report.tsv")),
                              separateMassAcc,
                              matrices,
                              carbamidomethyl,
                              noProtInf,
                              relaxProtInf,
                              peakCenter,
                              massAcc,
                              reportLibInfo))

}
