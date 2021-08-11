#' Create a predicted spectral library
#' @param msLib input MSLibrarian object
#' @param bestCe absolute path to a RData object containing optimal collision energies for MS/MS predictions.
#' @param predDb character giving the absolute location of the prediction database.
#' @param ceMode character specifying the mode for collision energy optimization. Available modes are c("length_only", "length_charge", "charge", scalar)
#' @export make.spectra.lib

make.spectra.lib <- function(msLib, bestCe, predDb, ceMode) {

  get.pept.length <- function(selCe, ce) {

    ce$peptideLength[ce$ce == selCe]


  }
  set.charge <- function(charge, bestCe) {

    ce = bestCe[bestCe$charge == charge, ]
    peptLengths = sapply(unique(ce$ce), get.pept.length, ce)
    if(length(unique(ce$ce)) == 1) {
      peptLengths = list(peptLengths)
    }
    names(peptLengths) = paste("ce_",unique(ce$ce), sep = "")
    peptLengths

  }

  get.msms.idx <- function(x) {

    seq(x["startIdx"], x["endIdx"])

  }
  pept.length.idx = function(aaLength,dbPeptLength) {

    which(dbPeptLength == aaLength)

  }
  ceTmp = lapply(unique(bestCe$charge), set.charge, bestCe)
  names(ceTmp) = paste("charge_",unique(bestCe$charge), sep = "")
  bestCe = ceTmp
  rm(ceTmp)
  gc()

  # ceMode = c("length_only", "length_charge", "charge", "fixed")[1]
  # fixedVal = 32
  msLib@PredLib@PrositLib = list()
  print("Connecting to Prosit database")
  # prositdb = dbConnect(RSQLite::SQLite(), "D:/Data_PROSIT/Libraries/Mouse_Yeast/SQLITE/mus_sac_comb_prositdb_2020.sqlite")
  prositdb = dbConnect(RSQLite::SQLite(), predDb)
  print("Importing precursor data from Prosit database...")
  precursorData = dbGetQuery(conn = prositdb, "SELECT * FROM PrecursorData")
  print("Adding matched precursors to slot PredLib@PrositLib...")
  msLib@PredLib@PrositLib$PrecursorData = precursorData[msLib@Sequences@Precursors$MzPass$predIdx,]

  if(identical(msLib@Sequences@Precursors$MzPass$peptide_sequence, msLib@PredLib@PrositLib$PrecursorData$StrippedPeptide)) {
    print("All precursors match!")
    print("Adding Uniprot Ids...")
    msLib@PredLib@PrositLib$PrecursorData$UniprotId = msLib@Sequences@Precursors$MzPass$protein_id

  } else {
    print("Precursors do not match...terminating code execution.")
    stop()
  }
  if(ceMode == "length_only") {
    print("Peptide length optimized CE selected...")
    ce = bestCe$charge_length_only
  } else if(ceMode == "length_charge") {
    print("Peptide length and charge optimized CE selected...")
    ce = do.call('c', bestCe[!grepl("only", names(bestCe))])
  } else if (ceMode == "charge") {
    print("Precursor charge optimized CE selected")
    ce = do.call('c', bestCe[grep("charge_only", names(bestCe))])
  } else if(is.numeric(ceMode) & ceMode >= 20 & ceMode <= 40) {
    print(paste("Fixed CE =", as.integer(ceMode), "set for all precursors..."))
    ce = as.integer(ceMode)
    names(ce) = paste("ce_", as.integer(ceMode), sep = "")
  } else {
    print("ceMode does not exist...terminating code execution...")
    stop()
  }

  print("Importing MS/MS data from database...")
  precursorData$PeptideLength = sapply(precursorData$StrippedPeptide, nchar)
  relativeIntensity = vector("numeric", precursorData$endIdx[nrow(precursorData)])
  for (i in 1:length(ce)) {

    if(ceMode == "length_only") {
      print(paste("Get indices for precursors with", paste(ce[[i]], collapse = ", "), "AA in Prosit DB"))
      selIdx = unlist(apply(precursorData[unlist(lapply(ce[[i]], pept.length.idx, precursorData$PeptideLength)), c("startIdx", "endIdx")],1, get.msms.idx))
    } else if(ceMode == "length_charge") {
      print(paste("Get indices for precursors with", paste(ce[[i]], collapse = ", "), "AA and a charge of", str_extract(names(ce[i]), "\\d+"), "in Prosit DB"))
      selIdx = unlist(apply(precursorData[intersect(unlist(lapply(ce[[i]], pept.length.idx, precursorData$PeptideLength)),
                                                    which(precursorData$PrecursorCharge == as.numeric(str_extract(names(ce[i]), "\\d+")))),
                                          c("startIdx", "endIdx")],1, get.msms.idx))

    } else if (ceMode == "charge") {
      print(paste("Getting all indices for CE =", str_remove(names(ce[i]), pattern = ".*ce_"), "in Prosit library"))
      selIdx = unlist(apply(precursorData[which(precursorData$PrecursorCharge == as.numeric(str_extract(names(ce[i]),
                                                                                                        "\\d+"))),
                                          c("startIdx", "endIdx")], 1, get.msms.idx))

    } else if (is.numeric(ceMode) & ceMode >= 20 & ceMode <= 40) {
      print(paste("Getting all indices for CE =", ce, "in Prosit library"))
      selIdx = unlist(apply(precursorData[,c("startIdx", "endIdx")], 1, get.msms.idx))

    }

    print(paste("Selecting relative intensities predicted with a CE =", str_remove(names(ce[i]), pattern = ".*ce_")))
    relativeIntensity[selIdx]= dbGetQuery(conn = prositdb, paste("SELECT RelativeIntensity_",str_remove(names(ce[i]), pattern = ".*ce_"), " FROM MsmsData", sep = ""))[selIdx,]
    print("Subprocess completed!")

  }
  rm(selIdx, precursorData)
  gc()
  print("Getting MS/MS indices for matched precursors...")
  msmsIdx = apply(msLib@PredLib@PrositLib$PrecursorData[, c("startIdx", "endIdx")], 1, get.msms.idx)
  endIdx = cumsum(unlist(lapply(msmsIdx, length)))
  msLib@PredLib@PrositLib$PrecursorData$startIdx = c(1, endIdx[1:(length(endIdx)-1)]+1)
  msLib@PredLib@PrositLib$PrecursorData$endIdx = endIdx
  rm(endIdx)

  print("Adding MS/MS data for matched precursors to slot PredLib@PrositLib...")
  msLib@PredLib@PrositLib$MsmsData = dbGetQuery(conn = prositdb, paste("SELECT FragmentMz, FragmentNumber, FragmentType, FragmentCharge FROM MsmsData"))
  msLib@PredLib@PrositLib$MsmsData$RelativeIntensity = relativeIntensity
  rm(relativeIntensity)
  msLib@PredLib@PrositLib$MsmsData = msLib@PredLib@PrositLib$MsmsData[unlist(msmsIdx),]
  rownames(msLib@PredLib@PrositLib$MsmsData) = 1:nrow(msLib@PredLib@PrositLib$MsmsData)
  print("Removing intensities equal to zero...")
  noZero = setdiff(which(msLib@PredLib@PrositLib$MsmsData$RelativeIntensity == 0),
                   which(msLib@PredLib@PrositLib$MsmsData$FragmentType == "y" &
                           msLib@PredLib@PrositLib$MsmsData$FragmentNumber == 1 &
                           msLib@PredLib@PrositLib$MsmsData$FragmentCharge == 1))
  msLib@PredLib@PrositLib$MsmsData = msLib@PredLib@PrositLib$MsmsData[-noZero,]
  print("Re-indexing of MS/MS indices...")
  startIdx = which(msLib@PredLib@PrositLib$MsmsData$FragmentType == "y" &
                     msLib@PredLib@PrositLib$MsmsData$FragmentNumber == 1 &
                     msLib@PredLib@PrositLib$MsmsData$FragmentCharge == 1)
  msLib@PredLib@PrositLib$PrecursorData$endIdx = c(startIdx[2:length(startIdx)]-1, nrow(msLib@PredLib@PrositLib$MsmsData))
  msLib@PredLib@PrositLib$PrecursorData$startIdx = startIdx
  rm(startIdx, noZero, msmsIdx)
  gc()
  print("Processing completed!")
  dbDisconnect(prositdb)
  msLib
}
