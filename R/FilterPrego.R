#' Filter spectral library based on top N high responding peptides for each protein group
#' @param inputLib absolute path to input spectral library (Spectronaut .CSV format or OpenSwath format .TSV)
#' @param outputLib absolute path to the filtered output library
#' @param topPept integer giving the top N high responding peptides for each protein group in the library
#' @param pregoPath path to the classify.exe tool as part of the PREGO distribution
#' @param replace logical determining if the input library should be replaced by the filtered library  (default is FALSE)
#' @export filter.prego

filter.prego <- function(inputLib = NULL, outputLib = NULL, topPept = NULL, pregoPath = NULL, replace = F) {


  get.lib.idx <- function(mat) {
    seq(mat["startIdx"], mat["endIdx"])
  }
  if(is.null(topPept)) {
    stop("No filter set! Please set filter using argument top")
  }
  if(is.null(inputLib) | !file.exists(inputLib) | !str_detect(inputLib, ".tsv$|.csv$")) {
    stop("Argument 'inputLib' missing or file does not exist with the correct extensions (.tsv or .csv)")
  }
  if(is.null(pregoPath)) {
    print(str_c("Searching for PREGO executable (classify.exe)"))
    pregoPath = system2("where", args = c("/r", "C:\\", "classify.exe"), stdout = T)
    if(length(pregoPath) == 1) {
      print("Found PREGO executable: classify.exe...")
    } else {
      stop("Could not auto-detect the PREGO executable: classify.exe. Please add path to executable (pregoPath argument).")
    }
  } else if(!(file.exists(pregoPath) & str_detect(pregoPath, "classify.exe"))) {
    stop("Specified PREGO executable: classify.exe does not exist...")
  }
  if(is.null(outputLib)) {
    outputLib = file.path(str_replace(inputLib, ".tsv$|.csv$",str_c("_prego_top_",topPept,".", tools::file_ext(inputLib))))
  } else if(identical(inputLib, outputLib) & !replace) {
    stop("Overwriting the input library is forbidden when argument replace is set to false")
  }
  if(replace) {
    print("Arg - replace = TRUE")
    if(identical(inputLib, outputLib)) {
      print("Will overwrite the input library after filtering.")
      outputLib = inputLib
    } else {
      print("Will remove input library after processing...")
      remove = T
    }
  }
  if(!identical(file_ext(inputLib), file_ext(outputLib))) {
    stop("Output library must have the same extension as the input library")
  } else {
    if(!dir.exists(dirname(outputLib))) {
      stop("Output directory does not exist...")
    }
  }
  if(str_detect(inputLib, ".csv$")) {
    print(str_c("Importing spectronaut library : ", inputLib))
    specLib = read_csv(inputLib, col_types = cols_only(RelativeIntensity = col_double(),
                                                       ModifiedPeptide = col_character(),
                                                       PrecursorCharge = col_integer(),
                                                       StrippedPeptide = col_character(),
                                                       UniprotId = col_character()))
    delim = ","

  } else if(str_detect(inputLib, ".tsv$")) {
    print(str_c("Importing OpenSwath library: ", inputLib))
    specLib = read_tsv(inputLib, col_types = cols_only(LibraryIntensity = col_double(),
                                                       ModifiedPeptideSequence = col_character(),
                                                       PrecursorCharge = col_integer(),
                                                       PeptideSequence = col_character(),
                                                       UniprotId = col_character()))
    delim = "\t"
    colnames(specLib)[grep("LibraryIntensity", colnames(specLib))] = "RelativeIntensity"
    colnames(specLib)[grep("ModifiedPeptideSequence", colnames(specLib))] = "ModifiedPeptide"
    colnames(specLib)[grep("PeptideSequence", colnames(specLib))] = "StrippedPeptide"
  }
  print("Extracting precursors...")
  precursors = str_c(specLib$ModifiedPeptide, "_", specLib$PrecursorCharge)
  startIdx = which(!duplicated(precursors))
  print(str_c("Number of precursors in library: ", length(startIdx)))
  seqs = as.data.frame(specLib[startIdx, c("UniprotId", "StrippedPeptide", "PrecursorCharge")])
  seqs$startIdx = startIdx
  seqs$endIdx = c(startIdx[2:length(startIdx)] - 1, nrow(specLib))
  rm(specLib)
  gc()
  seqs$key = str_c(seqs$UniprotId, "_", seqs$StrippedPeptide, "_", seqs$PrecursorCharge)
  pregoList = list()
  print("Running PREGO to rank high responding peptides...")
  for (i in sort(unique(seqs$PrecursorCharge))) {

    pregoData = data.frame(names = seqs$UniprotId[seqs$PrecursorCharge == i],
                           sequences = seqs$StrippedPeptide[seqs$PrecursorCharge == i])
    pregoInput = str_c(str_extract(inputLib, ".*\\/"), "prego_input.txt")
    pregoOutput = str_replace(pregoInput, "_input", "_scores")

    write_delim(pregoData,
                file = pregoInput,
                col_names = F,
                delim = "\t",
                quote_escape = F)
    system2(pregoPath, args = c(pregoInput),
            invisible = F,
            minimized = F,
            stdout = pregoOutput)

    pregoOut = read_delim(pregoOutput, skip = 1, delim = "\t", col_names = F, col_types = cols())
    colnames(pregoOut) = c("Accession", "Rank", "Sequence", "Score")
    pregoOut$key = str_c(pregoOut$Accession, "_", pregoOut$Sequence, "_", i)
    pregoList[[str_c("charge_", i)]] = pregoOut
    rm(pregoOut)

  }
  print("Retrieving results...")
  pregoOut = as.data.frame(do.call("rbind", pregoList))
  rownames(pregoOut) = pregoOut$key
  pregoOut = pregoOut[seqs$key,]

  if(identical(seqs$key, pregoOut$key)) {
    seqs$rank = pregoOut$Rank
    print(str_c("Selecting top ", topPept, " high responding peptides for each protein group"))
    seqs = seqs[seqs$rank <= topPept,]
    print(str_c("Number of precursors in library after filtering: ", nrow(seqs)))
    libIdx = unlist(apply(seqs[,c("startIdx", "endIdx")], 1, get.lib.idx))
    print(str_c("Writing filtered library to: ", outputLib))
    write_delim(x = read_delim(inputLib, delim = delim, col_types = cols())[libIdx,], file = outputLib, delim = delim)
    gc()
  } else {
    stop("Sequences do not match...terminating code execution.")
  }
  if(remove) {
    file.remove(inputLib)
  }
  outputLib
}

