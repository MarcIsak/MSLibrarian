#' Returns peptide sequences resulting from protein digestion based on cleavage rules, max.missed cleavages,
#' @param msLib a MSLibrarian object
#' @param rowStr a string defining the row names.
#' @param enzyme the enzyme to be used for digestion. The enzyme cleaves proteins according to the specified cleavage rules.
#' @param maxMissed maximum missed clevages allowed for a peptide sequence
#' @param carbamidomethyl Determines if every Cysteine residue is carbamidomethylated(TRUE) or not (FALSE). If TRUE, every cysteine adds a mass of 57.021464.
#' @param threads number of threads to use for the computation
#' @export digest.proteins

digest.proteins = function(msLib, rowStr, enzyme,maxMissed,carbamidomethyl,threads){

  rep.acc <- function(acc) {

    rep(acc[1], acc[2])

  }

  get.masses =  function(seq, k) {

    peptideMass = sum(str_count(seq,as.character(get.amino.acids()$AA[2:nrow(get.amino.acids())]))*
                        get.amino.acids()$ResidueMass[2:nrow(get.amino.acids())]) +
      (2*get.atomic.mass()["H"] + get.atomic.mass()["O"]) + (str_count(seq,"C")*57.021464)*k

  }

  if(carbamidomethyl == T) {
    k = 1
    print("Assuming all Cysteines are Carbamidomethylated")
  } else {
    k = 0
  }

  print(paste("In-silico digestion of proteins into peptides using", enzyme))
  seqs = cleave(as.character(msLib@Sequences@Proteins$Sequences), enzym = enzyme, missedCleavages = maxMissed)
  names(seqs) = msLib@Sequences@Proteins[,rowStr]
  peptNbr = unlist(lapply(seqs, length))



  peptideData = data.frame(protein_id = unlist(apply(cbind(names(peptNbr), peptNbr), 1, rep.acc)),
                           peptide_sequence = unlist(seqs),
                           stringsAsFactors = F,
                           row.names = 1:length(unlist(seqs)))

  print("Calculating masses...")
  cl <- makeCluster(threads)
  clusterEvalQ(cl, {
    .libPaths(.libPaths())
    library(MSLibrarian)
  })
  peptideData$precursor_mass = unlist(parallel::parLapply(cl, peptideData$peptide_sequence, get.masses, k))
  stopCluster(cl)
  peptideData$peptide_length = sapply(peptideData$peptide_sequence, nchar)
  msLib@Sequences@Peptides$All = peptideData
  print("Done!")
  msLib
}
