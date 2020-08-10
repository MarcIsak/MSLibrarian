#' Generates precursor m/z values for a filtered list of peptides, given the charge range set.
#' @param msLib a scell object
#' @param mzRange sets the interval for the allowed m/z values.
#' @param chargeRange charge ranges for resulting precursor ions
#' @param threads number of threads to use for the computation
#' @param outFile absolute path to output file with console information from parallel workers
#' @export get.precursors.tmp

get.precursors.tmp = function(msLib, mzRange,chargeRange, threads, outFile) {

  names.rep = function(n, name) {
    counter <<- counter + 1
    rep(name[counter], n)
  }
  unique.seqs = function(nonRed, dupSeq){
    tmp = dupSeq[dupSeq[,1] == nonRed[1],2]
    paste(c(nonRed[2],tmp), collapse=";")
  }
  mzVals <- function(x,protonMass) {

    y = (x[1] + x[2]*protonMass)/x[2]
  }

  nPrec = unlist(lapply(msLib@Sequences@FilterPeptides$Sequences, length))
  counter <<- 0
  names = unlist(lapply(nPrec, names.rep, names(nPrec)))

  massVals = unlist(msLib@Sequences@FilterPeptides$Masses[!is.na(msLib@Sequences@FilterPeptides$Masses)])
  seqs = as.data.frame(cbind(unlist(msLib@Sequences@FilterPeptides$Sequences[!is.na(msLib@Sequences@FilterPeptides$Masses)]),names),
                       stringsAsFactors = F)
  colnames(seqs) = c("Sequence", "Names")
  massVals = massVals[!duplicated(seqs[,"Sequence"])]
  dupSeqs = seqs[duplicated(seqs[,"Sequence"]),]
  seqs = seqs[!duplicated(seqs[,"Sequence"]),]
  print("Removing duplicate peptide sequences")
  cl <- makeCluster(threads)
  names = parallel::parApply(cl, seqs, 1, unique.seqs, dupSeqs)
  stopCluster(cl)
  gc()

  protonMass = 1.007276499879 # Mass of a proton in Da

  tmp = cbind(rep(massVals, each = length(seq(min(chargeRange),
                                              max(chargeRange)))),
              rep(seq(min(chargeRange),max(chargeRange)),
                  times = length(seq(min(chargeRange),max(chargeRange)))))
  cl <- makeCluster(threads)
  print("Calculating M/Z values for each precursor")
  out = parallel::parApply(cl, tmp, 1, mzVals, protonMass) # Consider using in other functions as well. Really good!
  stopCluster(cl)

  msLib@Sequences@Precursors = as.data.frame(cbind(unlist(lapply(seqs[,"Sequence"],rep,length(seq(min(chargeRange),max(chargeRange))))),
                                      unlist(lapply(massVals,rep,length(seq(min(chargeRange),max(chargeRange))))),
                                      out,
                                      rep(seq(min(chargeRange),max(chargeRange)),nrow(seqs)),
                                      unlist(lapply(names,rep,length(seq(min(chargeRange),max(chargeRange)))))), stringsAsFactors = F)

  colnames(msLib@Sequences@Precursors) = c("Sequence", "Mass / Da", "M/Z", "Charge", "ID")
  print("Filtering the precursor list based on the allowed M/Z range...")
  passMassLim = as.numeric(msLib@Sequences@Precursors[,"M/Z"]) >= min(mzRange) &
    as.numeric(msLib@Sequences@Precursors[,"M/Z"]) <= max(mzRange)
  msLib@Sequences@Precursors = msLib@Sequences@Precursors[passMassLim,]
  print("Sorting the precursor list based on the charge state...")
  msLib@Sequences@Precursors = msLib@Sequences@Precursors[order(as.numeric(msLib@Sequences@Precursors[,"Charge"])),]
  print("Done!")
  msLib
  # list(seqs = unlist(lapply(seqs[,"Sequence"],rep,length(seq(min(chargeRange),max(chargeRange))))),
  #      massVals = unlist(lapply(massVals,rep,length(seq(min(chargeRange),max(chargeRange))))),
  #      mz = out,
  #      charges = rep(seq(min(chargeRange),max(chargeRange)),length(seqs)),
  #      names = unlist(lapply(names,rep,length(seq(min(chargeRange),max(chargeRange))))))

}


