#' Loads a processed Prosit Library file into a MSLibrarian object
#' @param msLib input MSLibrarian object
#' @param file absolute path to processed Prosit library file
#' @param chargeRange numeric vector containing the range of allowed charge states
#' @param mzRange numeric vector containing the range of allowed M/Z values
#' @export load.prosit.lib

load.prosit.lib <- function(msLib, file, chargeRange, mzRange) {

  load(file)

  predictSeqs = rep(msLib@Sequences@Peptides$Predictable$peptide_sequence, length(chargeRange))

  if(all(outLib$PrecursorData$StrippedPeptide == predictSeqs & length(outLib$PrecursorData$StrippedPeptide) == length(predictSeqs))) {

    print("Sequences match...")
    print(paste("Selecting top", as.character(max(msLib@Sequences@Peptides$TopPrego$prego_rank)), "PREGO hits..."))
    outLib$PrecursorData = outLib$PrecursorData[rep(msLib@Sequences@Peptides$Predictable$prego_pass, length(chargeRange)),]
    print("Removing duplicated sequences")
    outLib$PrecursorData = outLib$PrecursorData[!rep(msLib@Sequences@Peptides$TopPrego$duplicated_sequence, length(chargeRange)),]
    print(paste("Selecting precursors passing the M/Z filter:", as.character(min(mzRange)), "-", as.character(max(mzRange))))
    outLib$PrecursorData = outLib$PrecursorData[msLib@Sequences@Precursors$Unique$mz_pass,]
    if(all(outLib$PrecursorData$StrippedPeptide == msLib@Sequences@Precursors$MzPass$peptide_sequence &
       length(outLib$PrecursorData$StrippedPeptide) == length(msLib@Sequences@Precursors$MzPass$peptide_sequence))) {
      print("Precursors match...")
      print("Adding Uniprot IDs...")
      outLib$PrecursorData$UniprotId = msLib@Sequences@Precursors$MzPass$protein_id
      msLib@PredLib@PrositLib = outLib
      print("Done! Library added to slot PredLib@PrositLib.")
      msLib
    } else {
      print("Precursors do NOT match...terminating code execution.")
      rm(outLib)
      stop()
    }
  } else {
    print("Sequences do not match...terminating code execution.")
    rm(outLib)
    stop()
  }
}
