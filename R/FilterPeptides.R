#' Filters peptide sequences based on amino acid length
#' @param msLib a scell input object
#' @param AALengthRange minimum peptide length allowed
#' @param rejectAAs amino acids to reject. Peptides containing these amino acids will be removed.
#' @param mzRange mass over charge range to consider for resulting precursor ions
#' @param chargeRange charges allowed for resulting precursor ions.
#' @param threads number of threads to use for the computation
#' @param outFile absolute path to output file for parallel workers
#' @export filter.peptides

filter.peptides = function(msLib, AALengthRange, rejectAAs, mzRange, chargeRange, threads, outFile){

  get.seq <- function(x) {
    x$Seq
  }
  get.mass <- function(x) {
    x$Mass
  }
  protonMass = 1.007276499879 # Mass of a proton in Da
  massRange = c((min(mzRange)*min(chargeRange)) - min(chargeRange)*protonMass,
                 (max(mzRange)*max(chargeRange)) - max(chargeRange)*protonMass) # Allowed mass range for peptides

  #passList = vector("list",length(msLib@Sequences@FilterPeptides$Sequences))

  cl <- makeCluster(threads, outfile = outFile)
  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths("C:/Users/marc/Dropbox (Human Neural Develop)/Marc/MISTR/R/win-library/3.5"))

  filterList = foreach(i = 1:length(msLib@Sequences@Peptides$Sequences), .packages = "tcltk") %dopar% {

    if(!exists("pb")) {
      pb = tkProgressBar("Parallel task",min = 1, max = length(msLib@Sequences@Peptides$Sequences))
    }
    setTkProgressBar(pb, i)
    pass = apply(rbind(nchar(msLib@Sequences@Peptides$Sequences[[i]]) >= min(AALengthRange) & nchar(msLib@Sequences@Peptides$Sequences[[i]]) <= max(AALengthRange),
                       msLib@Sequences@Peptides$Masses[[i]] >= min(massRange) & msLib@Sequences@Peptides$Masses[[i]] <= max(massRange),
                       !apply(sapply(msLib@Sequences@Peptides$Sequences[[i]],str_count,rejectAAs),2,any)), 2, all)

    out = list(Seq = msLib@Sequences@Peptides$Sequences[[i]][pass],
               Mass = msLib@Sequences@Peptides$Masses[[i]][pass])
    #passList[[i]] = pass
  }
  stopCluster(cl)
  print("pass_1")
  msLib@Sequences@FilterPeptides$Sequences = lapply(filterList,get.seq)
  msLib@Sequences@FilterPeptides$Masses = lapply(filterList,get.mass)
  print("pass_2")
  names(msLib@Sequences@FilterPeptides) = names(msLib@Sequences@Peptides)
  names(msLib@Sequences@FilterPeptides$Masses) = names(msLib@Sequences@Peptides$Masses)
  names(msLib@Sequences@FilterPeptides$Sequences) = names(msLib@Sequences@Peptides$Sequences)
  print("pass_4")
  msLib@Sequences@FilterPeptides$Masses[lapply(msLib@Sequences@FilterPeptides$Masses,length) == 0] = NA # If a an entry has no peptides after filtering, it gets a NA value. Or we can not use View() to look at the list
  print("pass_5")
  #sco@Digestion = passList
  #print(paste("Number of Peptides: ",do.call(sum,lapply(sco@FilterPeptides$Sequences,length)),sep="")) # Add to the Digestion slot later - summarizes the digestion settings.
  msLib
}
