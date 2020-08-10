#' Generates precursor m/z values for a filtered list of peptides, given the charge range set.
#' @param msLib a scell object
#' @param mzRange sets the interval for the allowed m/z values.
#' @param chargeRange charge ranges for resulting precursor ions
#' @param threads number of threads to use for the computation
#' @param outFile absolute path to output file with console information from parallel workers
#' @export get.precursors

get.precursors = function(msLib, mzRange,chargeRange, threads, outFile) {

  mzVals <- function(x,protonMass) {

    y = (x[1] + x[2]*protonMass)/x[2]
  }
  massVals = unlist(msLib@Sequences@FilterPeptides$Masses[!is.na(msLib@Sequences@FilterPeptides$Masses)])
  seqs = unlist(msLib@Sequences@FilterPeptides$Sequences[!is.na(msLib@Sequences@FilterPeptides$Masses)])

  cl <- makeCluster(threads)
  print("Starting parallel computing again")
  protonMass = 1.007276499879 # Mass of a proton in Da

  tmp = cbind(rep(massVals, each = length(seq(min(chargeRange),
                                              max(chargeRange)))),
              rep(seq(min(chargeRange),max(chargeRange)),
                  times = length(seq(min(chargeRange),max(chargeRange)))))


  out = parallel::parApply(cl, tmp, 1, mzVals, protonMass) # Consider using in other functions as well. Really good!
  stopCluster(cl)
  msLib@Sequences@Precursors = cbind(unlist(lapply(seqs,rep,length(seq(min(chargeRange),max(chargeRange))))),
                                     unlist(lapply(massVals,rep,length(seq(min(chargeRange),max(chargeRange))))),
                                     out,rep(seq(min(chargeRange),max(chargeRange)),length(seqs)))
  colnames(msLib@Sequences@Precursors) = c("Sequence", "Mass / Da", "M/Z", "Charge")
  passMassLim = as.numeric(msLib@Sequences@Precursors[,"M/Z"]) >= min(mzRange) &
    as.numeric(msLib@Sequences@Precursors[,"M/Z"]) <= max(mzRange)
  msLib@Sequences@Precursors = msLib@Sequences@Precursors[passMassLim,]
  msLib

}


