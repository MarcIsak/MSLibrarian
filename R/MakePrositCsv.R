#' Generates a CSV file with precursor information, to be used for prediction of MS/MS spectra in Prosit.
#' @param fasta input fasta file
#' @param chargeRange Charge ranges considered for precursor ions.
#' @param ceRange Collision energy range for predictions.
#' @param prefix file prefix to use (like the species from which the fasta file is derived)
#' @param outputFolder full path to the output CSV file
#' @param threads number of threads to use
#' @export make.prosit.csv

make.prosit.csv <- function(fasta, chargeRange, ceRange, prefix, outputFolder, threads = detectCores()) {

  if(threads > detectCores()) {
    threads = detectCores()
  }
  if(threads > 8) {
    threads = 8
  } # Should be a thread check function since it is used in several places...
  print("Load FASTA...")
  msLib = read.fasta(fasta)
  Sys.sleep(2)
  msLib = digest.proteins(msLib = msLib,
                          rowStr = "Accession",
                          enzyme = "trypsin",
                          maxMissed = 0,
                          carbamidomethyl = T,
                          threads = threads)
  print("Filter peptides...")
  msLib = filter.peptides.tmp(msLib = msLib,
                              AALengthRange = c(7, 30),
                              rejectAAs = c("U", "O", "B", "J", "X", "Z"),
                              mzRange = c(0, Inf),
                              chargeRange = chargeRange)
  msLib = get.precursors(msLib = msLib,
                         mzRange = c(0, Inf),
                         chargeRange = chargeRange,
                         matchDb = F,
                         threads = threads)

  print("Preparing for writing files...")

  precursorSeq = msLib@Sequences@Precursors$MzPass$peptide_sequence
  precursorCharge = msLib@Sequences@Precursors$MzPass$precursor_charge

  ceSel = seq(min(ceRange), max(ceRange))

  for (i in ceSel) {

    prositInput = data.frame(modified_sequence = precursorSeq,
                             collision_energy = rep(i, length(precursorSeq)),
                             precursor_charge = precursorCharge)
    csvStr = str_c(prefix,"_digest_ce", as.character(i),".csv", sep = "")
    print(str_c("writing: ", csvStr))
    write_csv(prositInput, file = paste(outputFolder, csvStr, sep = ""))

  }
  write_delim(data.frame(pred_file = str_c(prefix,"_digest_ce", ceSel,".csv", sep = ""),
                         task_id = NA),
              file = str_c(outputFolder, "task_id.txt"),
              delim = "\t")
  metadata = file.path(outputFolder, str_c(prefix, "_metadata.RData"))
  print(str_c("Saving SQLite metadata to:", metadata))
  save(msLib,
       ascii = F,
       compress = T,
       file = metadata) # Should be saved automatically when running make.prosit.csv
}


