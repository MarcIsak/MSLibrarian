#' Runs PeptideProphet and iProphet on input Comet search engine results
#' @param interactPath Absolute path to the interactparser.exe
#' @param peptideProphetPath Absolute path to the peptideprophetparser.exe
#' @param interProphetPath Absolute path to the interprophetparser.exe
#' @param cometFiles input Comet database search result files (pepXML)
#' @param output Absolute path to the desired output folder
#' @param threads integer giving the number of threads to use
#' @export run.prophets

run.prophets <- function(interactPath, peptideProphetPath, interProphetPath, cometFiles, output, threads) {

  decoy = "DECOY=DECOY"
  threads = str_c("THREADS=", threads)
  # nonparam = "NONPARAM"
  interactFile = file.path(output, "interact.pep.xml")

  print("Merging Comet output files...")
  system2(interactPath, # str_c(tppPath, "InteractParser.exe")
          args = c(interactFile,
                   cometFiles),
          invisible = F,
          minimized = T,
          wait = T)
  print("Running PeptideProphet on the merged results file...")
  system2(peptideProphetPath, # str_c(tppPath, "PeptideProphetParser.exe")
          args = c(interactFile,
                   "PPM",
                   "RT",
                   decoy,
                   "NONPARAM"),
          invisible = F,
          minimized = T,
          wait = T)
  print("Running iProphet on the PeptideProphet results file...")
  system2(interProphetPath, # str_c(tppPath, "InterProphetParser.exe")
          args = c(threads,
                   "NOFPKM",
                   decoy,
                   interactFile,
                   str_replace(interactFile, "interact", "interact_ipro")),
          invisible = F,
          minimized = T,
          wait = T)

}



# tic()
# interactParserPath = "C:/TPP/bin/InteractParser.exe"
# interactFile = "D:/MS_files/MSLibrarian_DIA_test/interact.pep.xml"
# #cometFiles = "D:/MS_files/MSLibrarian_DIA_test/comet_files/*pep.xml"
# cometFiles = list.files("D:/MS_files/MSLibrarian_DIA_test/", pattern = "pep.xml", full.names = T)
#
# system2(interactParserPath, args = c(interactFile, cometFiles),invisible = F, minimized = T)
# toc()

# tic()
# peptideProphetPath = "C:/TPP/bin/PeptideProphetParser.exe"
# cometFiles = "D:/MS_files/MSLibrarian_DIA_test/comet_files/interact.pep.xml"
# rt = "RT" # Will perhaps reduce the number of deviating peptides used for calibration in DeepLC
# ppm = "PPM"
# decoy = "DECOY=DECOY"
# nonparam = "NONPARAM"
#
# system2(peptideProphetPath , args = c(cometFiles, ppm, rt, decoy, nonparam), invisible = F, minimized = T)
# toc()

# tic()
# interProphetPath = "C:/TPP/bin/InterProphetParser.exe"
# nThreads = 10
# threads = paste("THREADS=", nThreads, sep = "")
# noFPKM = "NOFPKM"
# decoy = "DECOY=DECOY"
# peptideProphetFile = "D:/MS_files/MSLibrarian_DIA_test/comet_files/interact.pep.xml"
# outFile = "D:/MS_files/MSLibrarian_DIA_test/comet_files/interact_ipro.pep.xml"
#
# system2(interProphetPath, args = c(threads, noFPKM, decoy, peptideProphetFile, outFile))
# toc()


