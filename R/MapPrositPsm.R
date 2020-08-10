#' Maps indices of Prosit entries to PSMs in a PeptideShaker object.
#' @param pepShkr A PeptideShaker oject
#' @param prositLib A processed Prosit library file
#' @export map.prosit.psm

map.prosit.psm <- function(pepShkr, prositLib) {


  pepShkr@psmsPass[["specLib"]] = pepShkr@psmsPass$PSMs[pepShkr@psmsPass$bestPsmRows,] # Non-redundant list of PSMs with the highest scores (i.e a spectral library)

  # Consider allocating memory for specLib when creating pepShkr object in PepShkrProcess
  counter <<- 0
  find.prosit.match = function(x,prositLib){

    counter <<- counter + 1
    print(counter)
    f = which(prositLib$Precursor$StrippedPeptide == x[1] &
                prositLib$Precursor$PrecursorCharge == x[2])[1] # This is a temporary thing to avoid getting a list out, remove [1] once the Isoleucine issue is resolved.

  }
  peptSeqs = pepShkr@idData@peptides@peptides$pepSeq[pepShkr@psmsPass$specLib$peptSeq_indx] # Consider adding this as column instead of peptSeq_indx
  pepShkr@psmsPass$specLib$prositIndx = as.numeric(apply(cbind(peptSeqs,pepShkr@psmsPass$specLib$chargeState),1,
                                                         find.prosit.match,prositLib)) # returns a vector with indices matching prositLib to non-redundant lite of PSMs with the highest scores.


  # pepShkr@psmsPass$specLib$prositIndx = as.numeric(apply(pepShkr@psmsPass$specLib[,c("peptide_ref","chargeState")],1,
  #                                                        find.prosit.match,prositLib)) # returns a vector with indices matching prositLib to non-redundant lite of PSMs with the highest scores.
  pepShkr

}
