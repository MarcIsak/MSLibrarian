#' Runs a DIA-NN analysis of the selected DIA-file and input library in TSV format
#' @param diannPath Absolute path to the DIA-NN executable
#' @param diaFile Absolute path to the input DIA raw file
#' @param libFile Absolute path to the input spectral library file (as given by the make.osw.lib function)
#' @param output Desired file name for DIA-NN reports (incl. absolute path)
#' @param threads Integer specifying the number of threads for the analysis
#' @param fasta Absolute path to a FASTA file for matching of peptides sequences to proteins
#' @param enzyme Character specifying the enzyme to use for in-silico digestion of sequences in the FASTA file
#' @param separateMassAcc should mass accuracy be determined for each individual MS file (logical)
#' @export run.diann

run.diann <- function(diannPath, diaFile, libFile, output = NULL, threads = detectCores(), fasta = NULL, enzyme, separateMassAcc = F) {

  if(separateMassAcc) {
    separateMassAcc = "--individual-mass-acc"
  } else {
    separateMassAcc = NULL
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
  if(str_detect(libFile, ".tsv$")) {
    # Should have some function that checks that the input library has the correct data...
    if(sum(str_detect(colnames(read_tsv(libFile, n_max = 1, col_types = cols())), "ModifiedPeptideSequence|LibraryIntensity|NormalizedRetentionTime")) == 3) {
      print("OpenSwath format selected")
      headers = paste(c("ModifiedPeptideSequence", "PrecursorCharge", "PrecursorMz", "NormalizedRetentionTime",
                        "ProductMz", "LibraryIntensity", "UniprotId", "*", "*", "*", "Decoy", "ProductCharge",
                        "FragmentType", "FragmentSeriesNumber"), collapse = ",")
    } else {
      stop("TSV-file is not an OpenSwath library!")
    }

  } else if(str_detect(libFile, ".csv$")) {
    # Should have some function that checks that the input library has the correct data...
    if(sum(str_detect(colnames(read_csv(libFile, n_max = 1, col_types = cols())), "ModifiedPeptide|RelativeIntensity|iRT")) == 3) {
      print("OpenSwath format selected")
      headers = paste(c("ModifiedPeptide", "PrecursorCharge", "PrecursorMz", "iRT",
                        "FragmentMz", "RelativeIntensity", "UniprotId", "*", "*", "*", "*", "FragmentCharge",
                        "FragmentType", "FragmentNumber"), collapse = ",")
    } else {
      stop("CSV-file is not a Spectronaut library!")
    }
    print("Spectronaut format selected")
  }
  print("Start analysis in DIA-NN!")
  system2(diannPath, args = c(paste("--f", diaFile),
                              paste("--lib", libFile),
                              paste("--threads", as.character(threads)),
                              fasta,
                              enzyme,
                              paste("--library-headers", headers),
                              paste("--out", file.path(output, "report.tsv")),
                              separateMassAcc,
                              noProtInf))

}
