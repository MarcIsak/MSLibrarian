#' Get indices for predicted PSMs that match PSMs in a SpectrastLibrary object
#' @param msLib input MSLibrarian object
#' @param threads number of processors to use for the parallel computing
#' @param outFile absolute path to the output file from the parallel workers
#' @export get.pred.idx

get.pred.idx <- function(msLib, threads, outFile) {

  pred.idx <- function(x) {

    which(msLib@PredLib@PrositLib$PrecursorData$LabeledPeptide == as.character(x[1]) &
            msLib@PredLib@PrositLib$PrecursorData$PrecursorCharge == as.numeric(x[2]))

  }
  slotStr = c("Consensus","Unfiltered")
  for (i in 1:length(slotStr)) {

    tmpPsms = slot(msLib@SpectrastLib,slotStr[i])@PrecursorData$FilterLib[,c("ModifiedPeptideSequence","PrecursorCharge")]

    cl <- makeCluster(threads, outfile = outFile)
    registerDoParallel(cl)
    clusterEvalQ(cl, .libPaths("C:/Users/marc/Dropbox (Human Neural Develop)/Marc/MISTR/R/win-library/3.5"))

    print(paste("Begin processing of ",slotStr[i]," Library!",sep=""))

    outSeq = foreach(j = 1:nrow(tmpPsms), .packages = "tcltk", .combine = 'c') %dopar% {

      if(!exists("pb")) {
        pb = tkProgressBar("Parallel task",min = 1, max = nrow(tmpPsms))
      }
      setTkProgressBar(pb, j)
      pred.idx(tmpPsms[j,])

    }
    stopCluster(cl)
    slot(msLib@SpectrastLib,slotStr[i])@PrecursorData$FilterLib$predIdx = outSeq
    print(paste("Processing of ",slotStr[i]," Library completed!",sep=""))

  }
  gc()
  msLib
}
