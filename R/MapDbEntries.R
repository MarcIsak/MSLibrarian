#' Map entries in the predicted spectral library database to the peptides in the MSLibrarian object
#' @param msLib input MSLibrarian object.
#' @param predDb character giving the absolute location of the prediction database.
#' @param chargeRange numeric vector giving the range of precursor charge states to consider.
#' @export map.db.entries

map.db.entries <- function(msLib, predDb, chargeRange) {

  # set.idx <- function(peptideSeq, dbPrecursors, charge) {
  #
  #   query = str_c(peptideSeq, "_", charge)
  #   dbPrecursors = dbPrecursors[query,]
  #   dbPrecursors$precursor = str_c(dbPrecursors$StrippedPeptide, "_", charge)
  #   if(identical(dbPrecursors$precursor[!is.na(dbPrecursors$idx)], query[!is.na(dbPrecursors$idx)])) {
  #     print(str_c("Found a DB match for ", nrow(dbPrecursors[!is.na(dbPrecursors$idx),]),
  #                 " queried precursors (out of ", length(query), ") with a charge = ", charge))
  #   } else {
  #     stop("No matching precursors in database...")
  #   }
  #   dbPrecursors$idx
  #
  # }
  set.idx <- function(peptideSeq, dbPrecursors, charge) {

    query = str_c(peptideSeq, "_", charge)
    dbPrecursors = dbPrecursors[query]
    if(identical(dbPrecursors$key, query)) {
      print(str_c("Found a DB match for ", nrow(dbPrecursors[!is.na(dbPrecursors$idx),]),
                  " queried precursors (out of ", length(query), ") with a charge = ", charge))
    } else {
      stop("No matching precursors in database...")
    }
    dbPrecursors$idx
  }
  print("Connecting to SQLite database with Prosit predictions...")
  prositdb = dbConnect(RSQLite::SQLite(), predDb)
  precursorData = dbGetQuery(conn = prositdb, "SELECT StrippedPeptide, PrecursorCharge FROM PrecursorData")
  dbDisconnect(prositdb)
  rownames(precursorData) = str_c(precursorData$StrippedPeptide, "_", precursorData$PrecursorCharge)
  precursorData$idx = seq(1,nrow(precursorData))
  precursorData$key = str_c(precursorData$StrippedPeptide, "_", precursorData$PrecursorCharge)
  precursorData = as.data.table(precursorData)
  setkey(precursorData, key)
  print("Mapping precursors from database...")

  msLib@Sequences@Peptides$Predictable$predIdx_z2 = set.idx(peptideSeq = msLib@Sequences@Peptides$Predictable$peptide_sequence,
                                                            dbPrecursors = precursorData,
                                                            charge = min(chargeRange)) # max(chargeRange) would have to be changed in more charge states are looked for...
  msLib@Sequences@Peptides$Predictable$predIdx_z3 = set.idx(peptideSeq = msLib@Sequences@Peptides$Predictable$peptide_sequence,
                                                            dbPrecursors = precursorData,
                                                            charge = max(chargeRange)) # max(chargeRange)  would have to be changed in more charge states are looked for...
  msLib
}
