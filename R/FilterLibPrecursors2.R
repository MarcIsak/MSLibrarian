#' Filters slot PrecursorData$Consensus, and returns the filtered Precursor information in PrecursorData$ConsensusFilter
#' @param calibrationLib input SpectrastLibrary object
#' @param exUniMod numeric vector with UniMod accession numbers. A Precursor having any of these modification will be excluded.
#' @param precAARange numeric vector specifying the allowed length range in AA, that a Precursor is allowed to have.
#' @param chargeRange allowed charge range for precursors.
#' @param exAAs character vector describing AAs that a Precursor cannot have in its sequence.
#' @export filter.lib.precursors2

filter.lib.precursors2 <- function(calibrationLib, exUniMod, precAARange, chargeRange,exAAs) {
  mat = calibrationLib@PrecursorData$Original

  print("Filtering precursors...")
  filter = cbind(unlist(lapply(str_locate_all(mat$ModifiedPeptideSequence,paste("\\(UniMod:",as.character(exUniMod),"\\)",sep="")),length)) == 0,
                 sapply(mat$PeptideSequence,str_count,"") >= min(precAARange) & sapply(mat$PeptideSequence,str_count,"") <= max(precAARange),
                 unlist(lapply(str_locate_all(mat$PeptideSequence,c("U","O")),length)) == 0,
                 mat$PrecursorCharge >= min(chargeRange) & mat$PrecursorCharge <= max(chargeRange))
  calibrationLib@PrecursorData$FilterLib = mat[apply(filter,1,all),]
  print(str_c("Number of precursors after filtering: ",
              nrow(calibrationLib@PrecursorData$FilterLib), " ( ",  round(nrow(calibrationLib@PrecursorData$FilterLib)*100/nrow(mat),2), "% )"))
  calibrationLib
}
