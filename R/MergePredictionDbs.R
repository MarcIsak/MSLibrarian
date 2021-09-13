#' Merges prediction databases allowing for creation of hybrid-species spectral libraries
#' @param inputDB character vector with paths to prediction databases that will be merged
#' @param outputDB output merged prediction database
#' @export merge.prediction.dbs

merge.prediction.dbs <- function(inputDB, outputDB) {

  get.msms.seq = function(idx) {
    seq(idx["startIdx"],idx["endIdx"])
  }
  print("hihi")

  for(i in inputDB) {
    print("Fetching precursor data...")
    tmpDb = dbConnect(RSQLite::SQLite(), i)
    tmpPrec = dbGetQuery(tmpDb, "SELECT * FROM PrecursorData")
    print(str_c("Number of precursors: ", nrow(tmpPrec)))
    print("Fetching MS/MS data...")
    tmpMsms = dbGetQuery(tmpDb, "SELECT * FROM MsmsData")
    dbDisconnect(tmpDb)
    print(str_c("Number of MS/MS entries: ", nrow(tmpMsms)))
    gc()
    if(i == inputDB[1]) {
      outDb <- dbConnect(RSQLite::SQLite(), outputDB)
      print(str_c("Adding precursor data to database: ", outputDB))
      RSQLite::dbWriteTable(conn = outDb, name = "PrecursorData", value = tmpPrec, append = F)
      rm(tmpPrec)
      print(str_c("Adding MS/MS data to database:", outputDB))
      RSQLite::dbWriteTable(conn = outDb, name = "MsmsData", value = tmpMsms, append = F)
      dbDisconnect(outDb)
      rm(tmpMsms)
      print("Subprocess completed!")
    } else {
      outDb <- dbConnect(RSQLite::SQLite(), outputDB)
      precData = dbGetQuery(outDb, "SELECT * FROM PrecursorData")
      precursors = str_c(precData$LabeledPeptide, "_", precData$PrecursorCharge)
      lastIdx = precData$endIdx[nrow(precData)]
      dbDisconnect(outDb)
      rownames(tmpPrec) = str_c(tmpPrec$LabeledPeptide, "_", tmpPrec$PrecursorCharge)
      tmpPrec$idx = seq(1,nrow(tmpPrec))
      rmIdx = tmpPrec[intersect(rownames(tmpPrec), precursors),"idx"]
      msmsIdx = unlist(apply(tmpPrec[rmIdx, c("startIdx", "endIdx")], 1, get.msms.seq))
      print(str_c("Number of duplicated entries to remove: ", length(rmIdx)))
      tmpPrec = tmpPrec[-rmIdx, ]
      if(length(intersect(rownames(tmpPrec), precursors)) == 0) {
        tmpMsms = tmpMsms[-msmsIdx,]
        gc()
        print(str_c("Number of precursors after filtering: ", nrow(tmpPrec)))
        print(str_c("Number of MS/MS entries after filtering: ", nrow(tmpMsms)))
        diff = cumsum(tmpPrec$endIdx + 1 - tmpPrec$startIdx) + 1
        print("Reindexing of MS/MS entries...")
        tmpPrec$startIdx = c(1, diff[1:length(diff)-1]) + lastIdx
        tmpPrec$endIdx = c(tmpPrec$startIdx[2:length(tmpPrec$startIdx)] - 1, nrow(tmpMsms) + lastIdx)
        tmpPrec$idx = NULL
        rownames(tmpPrec) = seq(1, nrow(tmpPrec))
        outDb <- dbConnect(RSQLite::SQLite(), outputDB)
        print("Adding precursor data to combined database...")
        RSQLite::dbWriteTable(conn = outDb, name = "PrecursorData", value = tmpPrec, append = T)
        print("Adding MS/MS data to combined database...")
        RSQLite::dbWriteTable(conn = outDb, name = "MsmsData", value = tmpMsms, append = T)
        dbDisconnect(outDb)
      } else {
        stop("Duplicated precursors cannot be removed.")
      }
    }
  }
}


