#' Returns peptide sequences resulting from protein digestion based on cleavage rules, max.missed cleavages,
#' @param msLib a MSLibrarian object
#' @param rowStr a string defining the row names.
#' @param enzyme the enzyme to be used for digestion. The enzyme cleaves proteins according to the specified cleavage rules.
#' @param maxMissed maximum missed clevages allowed for a peptide sequence
#' @param carbamidomethyl Determines if every Cysteine residue is carbamidomethylated(TRUE) or not (FALSE). If TRUE, every cysteine adds a mass of 57.021464.
#' @param threads number of threads to use for the computation
#' @param outFile absolute path to output file from parallel workers
#' @export digest.proteins

digest.proteins = function(sco,rowStr,enzyme,maxMissed,carbamidomethyl,threads, outFile){

  msLib@Sequences@Peptides = list(Sequences = cleave(as.character(msLib@Sequences@Proteins$Sequences),enzyme,maxMissed),
                                  Masses = vector("list",nrow(msLib@Sequences@Proteins)))
  names(msLib@Sequences@Peptides$Sequences) = as.character(msLib@Sequences@Proteins[,rowStr])
  print(paste("Number of Peptides: ",do.call(sum,lapply(msLib@Sequences@Peptides$Sequences,length)),sep="")) # Add to the Digestion slot later - summarizes the digestion settings.

  if(carbamidomethyl == T) {
    k = 1
  } else {
    k = 0
  }

  get.masses =  function(x) {

    peptideMass = sum(str_count(x,as.character(get.amino.acids()$AA[2:nrow(get.amino.acids())]))*
                        get.amino.acids()$ResidueMass[2:nrow(get.amino.acids())]) +
      (2*get.atomic.mass()["H"] + get.atomic.mass()["O"]) + (str_count(x,"C")*57.021464)*k


    names(peptideMass) = NULL
    peptideMass

  }
  cl <- makeCluster(threads, outfile = outFile)
  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths("C:/Users/marc/Dropbox (Human Neural Develop)/Marc/MISTR/R/win-library/3.5"))

  masses = foreach(i = 1:nrow(msLib@Sequences@Proteins), .packages = c("tcltk","MSnbase")) %dopar% {

    if(!exists("pb")) {
      pb = tkProgressBar("Parallel task",min = 1, max = nrow(msLib@Sequences@Proteins))
    }
    setTkProgressBar(pb, i)
    #get.masses(as.matrix(msLib@Sequences@Peptides$Sequences[[i]]))
    apply(as.matrix(msLib@Sequences@Peptides$Sequences[[i]]),1,get.masses)
    #sco@Peptides$Masses[[i]] = apply(as.matrix(sco@Peptides$Sequences[[i]]),1,get.masses) # Get the monoisotopic masses (in Daltons) for each peptide sequence of every protein and add to slot Peptides and list Masses.
    # This operation takes several minutes, consider rewriting.
  }
  msLib@Sequences@Peptides$Masses = masses
  names(msLib@Sequences@Peptides$Masses) = names(msLib@Sequences@Peptides$Sequences)
  stopCluster(cl)
  msLib
}
