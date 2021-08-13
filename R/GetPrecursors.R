#' Generates precursor m/z values for a filtered list of peptides, given the charge range set.
#' @param msLib input MSLibrarian object
#' @param mzRange sets the interval for the allowed m/z values.
#' @param chargeRange charge ranges for resulting precursor ions
#' @param matchDb logical indicating if matching to a prediction database will be performed
#' @param threads number of threads to use
#' @export get.precursors

get.precursors = function(msLib, mzRange,chargeRange, matchDb, threads) {

  get.protein.id = function(precursors, duplicatedSeq) {

    paste(c(precursors[1],
            duplicatedPeptides$protein_id[duplicatedPeptides$peptide_sequence == precursors[2]]),
          collapse = ";")

  }
  get.mz.vals = function(x, protonMass) {

    (x[1] + x[2]*protonMass)/x[2]

  }

  msLib@Sequences@Peptides$Predictable$duplicated_sequence =
    duplicated(msLib@Sequences@Peptides$Predictable$peptide_sequence)

  print("Finding duplicated peptides...")
  duplicatedPeptides = msLib@Sequences@Peptides$Predictable[msLib@Sequences@Peptides$Predictable$duplicated_sequence,]
  print("Creating precursor data...")
  print("Removing duplicated sequences...")
  #precursorCols = c("protein_id", "peptide_sequence", "precursor_mass")
  precursorCols = c("protein_id", "peptide_sequence", "precursor_mass", str_c("predIdx_z", chargeRange))
  precursors = msLib@Sequences@Peptides$Predictable[!msLib@Sequences@Peptides$Predictable$duplicated_sequence, precursorCols]
  if(matchDb) {
    predIdx = msLib@Sequences@Peptides$Predictable[!msLib@Sequences@Peptides$Predictable$duplicated_sequence, paste("predIdx_z", chargeRange, sep = "")]
  }
  # cl = makeCluster(threads)
  # clusterEvalQ(cl, {
  #   .libPaths(.libPaths()) # Sets the library path
  # })
  # precursors$protein_id = parallel::parApply(cl,
  #                                            precursors[,c("protein_id", "peptide_sequence")],
  #                                            1,
  #                                            get.protein.id,
  #                                            duplicatedPeptides)
  # stopCluster(cl)
  print("hihi")
  duplicatedPeptides = duplicatedPeptides[,c("predIdx_z2", "protein_id")]
  duplicatedPeptides = aggregate(duplicatedPeptides, by = list(idx = duplicatedPeptides$predIdx_z2), rbind)
  duplicatedPeptides$protein_id = do.call('rbind', lapply(duplicatedPeptides$protein_id, str_c, collapse = ";"))
  duplicatedPeptides = as.data.table(duplicatedPeptides)
  setkey(duplicatedPeptides, idx)
  duplicatedPeptides = duplicatedPeptides[.(precursors$predIdx_z2)]
  if(identical(duplicatedPeptides$idx, precursors$predIdx_z2)) {
    nas = !is.na(duplicatedPeptides$protein_id)
    precursors$protein_id[nas] = str_c(precursors$protein_id[nas], duplicatedPeptides$protein_id[nas], sep = ";")
  } else {
    stop("Prediction indices do not match...")
  }
  gc()
  precursors = do.call('rbind', replicate(length(chargeRange), precursors, simplify = F))
  precursors$precursor_charge = rep(chargeRange, each = (nrow(precursors)/length(chargeRange)))
  toc()
  print("Calculating m/z values...")
  protonMass = 1.007276499879 # Mass of a proton in Da
  precursors$precursor_mz = apply(precursors[,c("precursor_mass", "precursor_charge")], 1, get.mz.vals, protonMass)
  print(str_c("Number of unique precursors: ", nrow(precursors)))
  if(matchDb){
    print("Database matching enabled...")
    print("Adding prediction indices...")
    precursors$predIdx = as.vector(as.matrix(predIdx))
    predictable = !is.na(precursors$predIdx)
    print(str_c("Number of unique predictable precursors: ", sum(predictable), " (", round(sum(predictable)*100/nrow(precursors), 2), "% )"))
    precursors = precursors[predictable,]
  }
  precursors$mz_pass = precursors$precursor_mz >= min(mzRange) & precursors$precursor_mz <= max(mzRange)
  msLib@Sequences@Precursors$Unique = precursors
  msLib@Sequences@Precursors$MzPass = precursors[precursors$mz_pass,-which(grepl("mz_pass", colnames(precursors)))]
  print(str_c("Number of unique predictable precursors within M/Z range (", min(mzRange), " - ",max(mzRange),"): ", sum(precursors$mz_pass)))
  msLib

}


