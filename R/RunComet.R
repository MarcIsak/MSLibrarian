#' Runs the database search engine Comet on input spectral files
#' @param cometPath character giving the absolute path to the Comet executable
#' @param msFiles absolute path to the input MS files (mzXML)
#' @param fasta character giving the absolute path to the desired FASTA file
#' @param cometParams character giving the absolute path to the Comet parameters file
#' @param threads integer specifying the number of threads to use for the conversion
#' @export run.comet

run.comet <- function(cometPath, msFiles, fasta, cometParams, threads) {

  comet.cmd <- function(msFile, cometPath, cometParams, fasta) {

    system2(cometPath,
            args = c(fasta, cometParams, msFile),
            invisible = F,
            minimized = T,
            wait = T)

  }
  fasta = str_c("-D", fasta) # Use OpenMS FileMerger algorithm to combine fasta files
  cometParams = str_c("-P", cometParams)

  if(threads > length(msFiles)) {
    threads = length(msFiles)
    print(str_c("Setting threads to ",threads, "..."))
  }

  cl = makeCluster(threads)
  parallel::parLapply(cl,
                      msFiles,
                      comet.cmd,
                      cometPath,
                      cometParams,
                      fasta)
  stopCluster(cl)
  gc()

}



#
# tic()
# cometPath = "C:/TPP/bin/comet.exe"

# #fastaFile = paste("-D", "D:/Databases/uniprot-filtered-proteome_UP000000589+AND+reviewed_yes+AND+organism__Mus_can.fasta", sep = "")

# # pseudoMsMsFiles = "D:/MS_files/MSLibrarian_DIA_test/comet_files/*.mzXML"
# pseudoMsMsFiles = list.files("D:/MS_files/MSLibrarian_DIA_test/comet_files/", pattern = "CK_",full.names = T)
# cl = makeCluster(6)
# parallel::parLapply(cl, pseudoMsMsFiles, run.comet, cometPath, cometParams, fastaFile)
# stopCluster(cl)
# gc()
# toc()
