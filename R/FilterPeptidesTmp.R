#' Filters peptide sequences based on amino acid length
#' @param msLib a scell input object
#' @param AALengthRange minimum peptide length allowed
#' @param rejectAAs amino acids to reject. Peptides containing these amino acids will be removed.
#' @param mzRange mass over charge range to consider for resulting precursor ions
#' @param chargeRange charges allowed for resulting precursor ions.
#' @export filter.peptides.tmp

filter.peptides.tmp = function(msLib, AALengthRange, rejectAAs, mzRange, chargeRange){

  protonMass = 1.007276499879 # Mass of a proton in Da
  massRange = c((min(mzRange)*min(chargeRange)) - min(chargeRange)*protonMass,
                (max(mzRange)*max(chargeRange)) - max(chargeRange)*protonMass) # Allowed mass range for peptides

  pass = cbind(msLib@Sequences@Peptides$All$precursor_mass >= min(massRange) &
                 msLib@Sequences@Peptides$All$precursor_mass <= max(massRange),
               msLib@Sequences@Peptides$All$peptide_length >= min(AALengthRange) &
                 msLib@Sequences@Peptides$All$peptide_length <= max(AALengthRange),
               !grepl(paste(rejectAAs, collapse = "|"), msLib@Sequences@Peptides$All$peptide_sequence))
  msLib@Sequences@Peptides$All$filter_pass = apply(pass, 1, all)
  msLib@Sequences@Peptides$Predictable = msLib@Sequences@Peptides$All[msLib@Sequences@Peptides$All$filter_pass,
                                                            -which(colnames(msLib@Sequences@Peptides$All) == "filter_pass")]
  rownames(msLib@Sequences@Peptides$Predictable) = 1:nrow(msLib@Sequences@Peptides$Predictable)
  msLib

}
