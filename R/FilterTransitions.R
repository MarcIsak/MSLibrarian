#' Performs transition filtering of a library, either based on the top N most intense transitions, or with a defined relative intensity cutoff
#' @param inputLib path to library that will be filtered in Spectronaut (.csv) format.
#' @param outputLib path to transition filtered library (Spectronaut *.csv)
#' @param topTrans integer specifying the top n transitions to retain for each precursor in the library
#' @param cutoffTrans a transition relative intensity cutoff (range 0-1). A cutoff = 0.05 means that transitions having an intensity < 5 % of the maximum intensity (100%) will be removed.
#' @param minTrans integer giving the minimum number of transitions that a precursor must have. Default is NULL, meaning that no filtering is performed.
#' @param replace logical indicating if the input library be replaced (default = FALSE)
#' @param threads number of threads to use for processing
#' @export filter.transitions

filter.transitions <- function(inputLib = NULL, outputLib = NULL, topTrans = NULL, cutoffTrans = NULL, minTrans = NULL, replace = F, threads = detectCores()) {

  sort.transitions <- function(mat, relativeIntensity) {

    order(relativeIntensity[seq(mat["startIdx"],
                                mat["endIdx"])], decreasing = T) + (mat["startIdx"] - 1)

  }
  top.n.transitions <- function(element,n) {

    element[seq(1,n)]

  }
  get.seq <- function(mat) {
    seq(mat["startIdx"], mat["endIdx"])
  }
  if(all(is.null(topTrans), is.null(cutoffTrans), is.null(minTrans))) {
    stop("No filter set! Please set filters using arguments topTrans, cutoffTrans and minTrans")
  }
  if(is.null(inputLib) | !file.exists(inputLib) | !str_detect(inputLib, ".tsv$|.csv$")) {
    stop("Argument 'inputLib' missing with the correct extensions (.tsv or .csv)")
  }
  defaultName  = F
  remove = F
  if(is.null(outputLib)) {
    outputLib = outputLib = file.path(str_replace(inputLib, ".tsv$|.csv$",str_c("_topN_", topTrans, "_cutoff_", cutoffTrans,"_trans_", minTrans, ".", tools::file_ext(inputLib))))
    defaultName = T
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
                                                       PrecursorCharge = col_integer()))
    delim = ","

  } else if(str_detect(inputLib, ".tsv$")) {
    print(str_c("Importing OpenSwath library: ", inputLib))
    specLib = read_tsv(inputLib, col_types = cols_only(LibraryIntensity = col_double(),
                                                       ModifiedPeptideSequence = col_character(),
                                                       PrecursorCharge = col_integer()))
    delim = "\t"
    colnames(specLib)[grep("LibraryIntensity", colnames(specLib))] = "RelativeIntensity"
    colnames(specLib)[grep("ModifiedPeptideSequence", colnames(specLib))] = "ModifiedPeptide"
  }
  precursors = str_c(specLib$ModifiedPeptide, "_", specLib$PrecursorCharge)
  startIdx = which(!duplicated(precursors))
  print(str_c("Found ", length(startIdx), " precursors and ", length(precursors)," transitions in the input library..."))
  relativeIntensity = specLib$RelativeIntensity
  rm(specLib)
  gc()
  endIdx = as.integer(c(startIdx[2:length(startIdx)] - 1, length(relativeIntensity)))
  pass = NULL
  if(!is.null(topTrans)) {
    topNValid = topTrans == floor(topTrans) & topTrans == ceiling(topTrans) & topTrans >= 1
    if(topNValid) {
      print(str_c("Selecting top ", topTrans, " most intense transitions for each precursor..."))
      cl = parallel::makeCluster(threads)
      clusterEvalQ(cl, {
        .libPaths(.libPaths())
      })
      transitions = parApply(cl, cbind(startIdx, endIdx), 1, sort.transitions, relativeIntensity)
      stopCluster(cl)
      idx = which(unlist(lapply(transitions, length)) > topTrans)
      if(length(idx) >= 1) {
        print(str_c("Performing transition filtering on: ", length(idx), " entries ( ", round(length(idx)*100/length(transitions),2), " % )"))
        transitions[idx] = lapply(transitions[idx], top.n.transitions, topTrans)
      } else {
        print(str_c("No precursor has more than ", topTrans, " transitions. No filtering will be performed..."))
      }
      pass = sort(unique(unlist(transitions))) # Does not require y1^1 ions to be present. New addition
      print(str_c("Library now contains: ", length(pass), " transitions ( ", round(length(pass)*100/length(relativeIntensity),2), "% of initial library )"))
      rm(transitions, idx)
    } else {
      stop("Invalid argument - topN")
    }
  } else {
    print("No topN filter set - skipping...")
    if(defaultName) {
      outputLib = str_remove(outputLib, "_topN_")
    }
  }
  if(!is.null(cutoffTrans)) {
    cutoffValid = is.numeric(cutoffTrans) & !is.integer(cutoffTrans) & between(cutoffTrans, 0,1)
    if(cutoffValid) {
      print(str_c("Removes transitions with a relative intensity < ", cutoffTrans))
      cutoffPass = unique(c(which(relativeIntensity >= cutoffTrans))) # Does not require the y1^1 ions to be present.
      if(is.null(pass)) {
        pass = cutoffPass
      } else {
        pass = sort(intersect(pass, cutoffPass))
      }
      rm(cutoffPass)
      print(str_c("Filtered library contains: ", length(pass), " transitions ( ", round(length(pass)*100/length(relativeIntensity),2), "% of initial library )"))
    } else {
      stop("Invalid argument - cutoff")
    }
  } else {
    print("No cutoff filter set - skipping...")
    if(defaultName) {
      outputLib = str_remove(outputLib, "_cutoff_")
    }
  }
  if(!is.null(minTrans)) {
    minTransValid = minTrans == floor(minTrans) & minTrans == ceiling(minTrans) & minTrans > 0
    if(minTransValid) {
      print(str_c("Removing all transitions for precursors having less than ", minTrans, " transitions."))
      if(!is.null(pass)) {
        precursors = precursors[pass]
        startIdx = which(!duplicated(precursors))
        endIdx  = c(startIdx[2:length(startIdx)] - 1, length(precursors))
      }
        diff = which((endIdx + 1 - startIdx) >= minTrans)
        minTransPass = unlist(apply(cbind(startIdx, endIdx)[diff,], 1, get.seq))
        print(str_c("Library now contains: ", length(minTransPass), " transitions ( ", round(length(minTransPass)*100/length(relativeIntensity),2), "% of initial library )"))
    } else {
      stop("Invalid argument - minTrans.")
    }
  } else {
    print("No minTrans filter set - skipping...")
    if(defaultName) {
      outputLib = str_remove(outputLib, "_trans_")
    }
  }
  rm(startIdx, endIdx, relativeIntensity, precursors)
  gc()
  specLib = read_delim(inputLib, delim = delim, col_types = cols())
  if(!is.null(pass)) {
    specLib = specLib[pass,]
    gc()
  }
  if(!is.null(minTrans)) {
    specLib = specLib[minTransPass,]
    rm(minTransPass)
  }
  print(str_c("Writing library to: ", outputLib))
  write_delim(specLib, file = outputLib, delim = delim)
  rm(specLib, pass)
  gc()
  if(remove) {
    file.remove(inputLib)
  }
  outputLib
}
