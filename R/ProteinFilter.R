#' Subsets an input spectral library based on protein IDs
#' @param inputLib path to input spectral library to be filtered
#' @param outputLib desired path to subsetted spectral library
#' @param proteinIds character vector of protein IDs which the subsetted library must contain
#' @export protein.filter

protein.filter <- function(inputLib, outputLib, proteinIds) {

  rep.idx <- function(mat) {
    rep(mat["rowId"], mat["nbrOfIds"])
  }
  as.seq <- function(mat) {
    seq(mat["startIdx"], mat["endIdx"])
  }
  if(str_detect(inputLib, ".csv$")) {
    print(str_c("Importing spectronaut library : ", inputLib))
    specLib = read_csv(inputLib, col_types = cols_only(UniprotId = col_character(),
                                                       ModifiedPeptide = col_character(),
                                                       PrecursorCharge = col_integer(),
                                                       StrippedPeptide = col_character()))
    delim = ","

  } else if(str_detect(inputLib, ".tsv$")) {
    print(str_c("Importing OpenSwath library: ", inputLib))
    specLib = read_tsv(inputLib, col_types = cols_only(UniprotId = col_character(),
                                                       ModifiedPeptideSequence = col_character(),
                                                       PrecursorCharge = col_integer(),
                                                       PeptideSequence = col_character()))
    delim = "\t"
    colnames(specLib)[grep("ModifiedPeptideSequence", colnames(specLib))] = "ModifiedPeptide"
    colnames(specLib)[grep("PeptideSequence", colnames(specLib))] = "StrippedPeptide"
  }
  startIdx = which(!duplicated(str_c(specLib$ModifiedPeptide, "_", specLib$PrecursorCharge)))
  endIdx = c(startIdx[2:length(startIdx)] -1, nrow(specLib))
  specLib = specLib[startIdx,]
  specLib$startIdx = startIdx
  specLib$endIdx = endIdx
  specLib[,c("PrecursorCharge", "StrippedPeptide", "ModifiedPeptide")] = NULL
  gc()
  specLib$rowId = seq_len(nrow(specLib))
  specLib$nbrOfIds = unlist(lapply(str_split(specLib$UniprotId, ";"), length))
  print(str_c("Unique protein groups in input spectral library: ", length(unique(specLib$UniprotId))))

  protTable = as.data.table(cbind(unlist(str_split(specLib$UniprotId, ";")), unlist(apply(specLib[,c("rowId", "nbrOfIds")], 1, rep.idx))))
  colnames(protTable) = c("UniprotId", "idx")
  setkey(protTable, UniprotId)

  #print(str_c("Unique protein groups in Calibration library: ", length(unique(calibLib@RawLib$ProteinId))))
  #proteinIds = unique(unlist(str_split(calibLib@RawLib$ProteinId, ";"))) # This data could be imported
  #proteinIds = str_remove(str_remove(proteinIds, "sp\\||tr\\|"), "\\|.*")
  #print(str_c("Unique protein Ids in Calibration library: ", length(proteinIds)))
  protTable = protTable[proteinIds]
  matchRow = sort(unique(as.numeric(protTable[!is.na(idx), idx])))
  acc = specLib[matchRow,"UniprotId"]
  print(str_c("Number of matching protein groups in the input library: ",nrow(unique(acc)), " ( ", round(nrow(unique(acc))*100/length(unique(specLib$UniprotId)),2), "% )"))
  print("Subsetting library")
  idx = unlist(apply(specLib[matchRow,c("startIdx", "endIdx")], 1, as.seq))
  write_delim(read_delim(inputLib, delim, col_types = cols())[idx,], file = outputLib, delim = delim)
  print(str_c("Library written to: ", outputLib))

}
