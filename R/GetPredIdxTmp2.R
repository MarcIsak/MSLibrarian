#' Get indices for predicted PSMs that match PSMs in a SpectrastLibrary object
#' @param msLib input MSLibrarian object
#' @param predStr character specifying the selected predicted library to find indices for
#' @param threads number of processors to use for the parallel computing
#' @param outFile absolute path to the output file from the parallel workers
#' @export get.pred.idx.tmp2

get.pred.idx.tmp2 <- function(msLib, predStr, threads, outFile) {

  pred.idx <- function(x, predStr) {

    tmp = which(slot(msLib@PredLib, selPred)$PrecursorData[,matchCol] == as.character(x[1]) &
            slot(msLib@PredLib, selPred)$PrecursorData$PrecursorCharge == as.numeric(x[2]))
    if (length(tmp) == 0) {
      tmp = NA
    }
    tmp
  }
  if(predStr == "prosit") {
    selPred = "PrositLib"
    predCol = "prositIdx"
    matchCol = "LabeledPeptide"
  } else if (predStr == "ms2pip") {
    selPred = "Ms2pipLib"
    predCol = "ms2pipIdx"
    matchCol = "StrippedPeptide"
  } else {
    print("Incorrect predicted library selected, check for typos!")
  }

  #expPsms = msLib@SpectrastLib@Unfiltered@PrecursorData$FilterLib[,c("ModifiedPeptideSequence","PrecursorCharge")]
  # idxCol = rep(NA, nrow(expPsms))
  cl <- makeCluster(threads, outfile = outFile)
  # registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths("C:/Users/marc/Dropbox (Human Neural Develop)/Marc/MISTR/R/win-library/3.5"))

  print(paste("Begin processing of the Unfiltered Library!",sep=""))

  outSeq = parallel::parApply(cl,
                              msLib@SpectrastLib@Unfiltered@PrecursorData$FilterLib[,c("ModifiedPeptideSequence","PrecursorCharge")],
                              1, pred.idx)

  # outSeq = foreach(j = 1:nrow(expPsms), .packages = "tcltk", .combine = 'c') %dopar% {
  #
  #   if(!exists("pb")) {
  #     pb = tkProgressBar("Parallel task",min = 1, max = nrow(expPsms))
  #   }
  #   setTkProgressBar(pb, j)
  #   tmp = pred.idx(expPsms[j,])
  #
  #   idxCol[j] = pred.idx(expPsms[j,])
  #
  # }
  stopCluster(cl)
  msLib@SpectrastLib@Unfiltered@PrecursorData$FilterLib[predCol] = outSeq
  print(paste("Processing of Unfiltered Library completed!",sep=""))
  gc()
  msLib

}
