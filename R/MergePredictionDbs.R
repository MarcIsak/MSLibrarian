#' Merges prediction databases allowing for creation of hybrid-species spectral libraries
#' @param inputDB character vector with paths to prediction databases that will be merged
#' @param outputDB output merged prediction database
#' @export merge.prediction.dbs

merge.prediction.dbs <- function(inputDB, outputDB) {

  get.prec.idx = function(seq, tmpDb) {
    which(tmpDb$StrippedPeptide == seq)

  }
  get.msms.seq = function(idx) {
    seq(idx["startIdx"],idx["endIdx"])
  }

  # inputDB = c("D:/Data_PROSIT/Libraries/Mouse/SQLITE/mus_musculus_prositdb_2020.sqlite",
  #             "D:/Data_PROSIT/Libraries/Yeast/SQLITE/saccharomyces_cerevisiae_prositdb_2020.sqlite")
  #
  # outputDb = "D:/Data_PROSIT/Libraries/Mouse_Yeast/SQLITE/mus_sac_comb_prositdb_2020.sqlite"

  for (i in inputDB) {
    print("Fetching precursor data...")
    tmpDb <- dbConnect(RSQLite::SQLite(), i)
    tmpPrec = dbGetQuery(tmpDb, "SELECT * FROM PrecursorData")
    print(str_c("Number of precursors: ", nrow(tmpPrec)))
    print("Fetching MS/MS data...")
    tmpMsms = dbGetQuery(tmpDb, "SELECT * FROM MsmsData")
    print(str_c("Number of MS/MS entries: ", nrow(tmpMsms)))
    gc()
    if(i == min(inputDB)) {
      outDb <- dbConnect(RSQLite::SQLite(), outputDB)
      print("Adding precursor data to combined database...")
      RSQLite::dbWriteTable(conn = outDb, name = "PrecursorData", value = tmpPrec, append = F)
      uniqueSeq = unique(tmpPrec$StrippedPeptide)
      precRows = nrow(tmpPrec)
      rm(tmpPrec)
      msmsRows = nrow(tmpMsms)
      print("Adding MS/MS data to combined database...")
      RSQLite::dbWriteTable(conn = outDb, name = "MsmsData", value = tmpMsms, append = F)
      dbDisconnect(tmpDb)
      print("Subprocess completed!")
    } else {
      rmIdx = sort(as.vector(sapply(intersect(uniqueSeq, tmpPrec$StrippedPeptide),
                                    get.prec.idx,
                                    tmpPrec)))
      print(str_c("Number of duplicated entries found: ", length(rmIdx)))
      msmsIdx = unlist(apply(tmpPrec[rmIdx, c("startIdx", "endIdx")], 1, get.msms.seq))
      print("Removing duplicated entries...")
      tmpPrec = tmpPrec[-rmIdx,]
      tmpMsms = tmpMsms[-msmsIdx, ]
      print(str_c("Number of precursors after filtering: ", nrow(tmpPrec)))
      print(str_c("Number of MS/MS entries after filtering: ", nrow(tmpMsms)))
      gc()
      rownames(tmpPrec) = seq(1, nrow(tmpPrec)) + precRows
      rownames(tmpMsms) = seq(1, nrow(tmpMsms)) + msmsRows
      print("Reindexing of MS/MS entries...")
      startIdx = which(tmpMsms$FragmentType == "y" & tmpMsms$FragmentNumber == 1 & tmpMsms$FragmentCharge == 1)
      tmpPrec$endIdx = c(startIdx[2:length(startIdx)]-1, nrow(tmpMsms)) + msmsRows
      tmpPrec$startIdx = startIdx + msmsRows
      print("Adding precursor data to combined database...")
      RSQLite::dbWriteTable(conn = outDb, name = "PrecursorData", value = tmpPrec, append = T)
      print("Adding MS/MS data to combined database...")
      RSQLite::dbWriteTable(conn = outDb, name = "MsmsData", value = tmpMsms, append = T)
      dbDisconnect(tmpDb)
    }
  }
  dbDisconnect(outDb)

}


