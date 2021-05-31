#' Creates a generic spectral library from predicted data, which can be used in OpenSWATH (TSV-format)
#' @param msLib input MSLibrarian object
#' @param format character defining the output library format c("spectronaut", "openswath")
#' @param calibrationLib path to calibration library
#' @param outputLib character specifying the output library file name (incl. absolute path)
#' @param batchSize numeric setting the maximum number of entries to process at a time (to keep the memory usage low)
#' @export export.spectra.lib

export.spectra.lib <- function(msLib, format, calibrationLib, outputLib, batchSize) {

  replicate.precursor.data <- function(precursorData, col) {

    rep(precursorData[col], as.numeric(precursorData["endIdx"]) - as.numeric(precursorData["startIdx"]) + 1)

  }

  get.idx <- function(precursorData) {

    seq(precursorData["startIdx"], precursorData["endIdx"])

  }

  predLibCols = c("PrecursorMz", "PrecursorCharge", "iRT", "StrippedPeptide", "LabeledPeptide", "UniprotId")

  if(format == "spectronaut") {

    libCols = c("RelativeIntensity", "FragmentMz", "ModifiedPeptide", "LabeledPeptide", "StrippedPeptide", "PrecursorCharge",
                        "PrecursorMz", "iRT", "FragmentNumber", "FragmentType", "FragmentCharge", "FragmentLossType", "UniprotId")
    selPrecCols = libCols[c(7, 6, 8, 5, 4, 13)]
    predLibCols = predLibCols[1:length(selPrecCols)]
    print("Spectronaut format selected...")
    if(!str_detect(outputLib, ".csv$")) {
      outputLib = str_c(outputLib, ".csv")
    }

  } else if(format == "openswath") {

    libCols = c("PrecursorMz", "ProductMz", "PrecursorCharge", "ProductCharge", "LibraryIntensity", "NormalizedRetentionTime",
                   "PeptideSequence", "ModifiedPeptideSequence", "PeptideGroupLabel", "LabelType", "CompoundName", "SumFormula",
                   "SMILES", "ProteinId", "UniprotId", "FragmentType", "FragmentSeriesNumber", "Annotation", "CollisionEnergy",
                   "PrecursorIonMobility", "TransitionGroupId", "TransitionId", "Decoy", "DetectingTransition", "IdentifyingTransition",
                   "QuantifyingTransition", "Peptidoforms")
    selPrecCols = libCols[c(1,3,6,7,8,15)]
    predLibCols = predLibCols[1:length(selPrecCols)]
    print("OpenSwath format selected...")
    if(!str_detect(outputLib, ".tsv$")) {
      outputLib = str_c(outputLib, ".tsv")
    }

  } else {
    print("Incorrect format specified...terminating code execution")
    stop()
  }

  if(nrow(msLib@PredLib@PrositLib$PrecursorData) <= batchSize) {

    batchSize = nrow(msLib@PredLib@PrositLib$PrecursorData)
    print(paste("Setting batchSize to",nrow(msLib@PredLib@PrositLib$PrecursorData)))
    nbrOfAppends = 1
    firstIdx = 1
    lastIdx = batchSize

  } else {
    nbrOfAppends = ceiling(nrow(msLib@PredLib@PrositLib$PrecursorData)/batchSize)
    firstIdx = c(1,seq(1, nbrOfAppends-1)*batchSize+1)
    lastIdx = c(seq(1,nbrOfAppends-1)*batchSize, nrow(msLib@PredLib@PrositLib$PrecursorData))
  }

  file.create(outputLib)

  for (i in 1:nbrOfAppends) {

    spectraLib = as.data.frame(matrix(NA,
                                  sum(msLib@PredLib@PrositLib$PrecursorData$endIdx[firstIdx[i]:lastIdx[i]] - msLib@PredLib@PrositLib$PrecursorData$startIdx[firstIdx[i]:lastIdx[i]] + 1),
                                  length(libCols)),
                           stringsAsFactors = F)
    colnames(spectraLib) = libCols
    gc()


    for (j in 1:length(selPrecCols)) {

      print(paste("Adding the", selPrecCols[j], "column to the", format, "library..."))
      spectraLib[, selPrecCols[j]] = unlist(apply(msLib@PredLib@PrositLib$PrecursorData[firstIdx[i]:lastIdx[i],], 1, replicate.precursor.data, predLibCols[j]))

    }
    gc()
    print(paste("Adding fragment ion data to the", format, "library..."))
    msmsIdx = unlist(apply(msLib@PredLib@PrositLib$PrecursorData[firstIdx[i]:lastIdx[i],], 1, get.idx))
    if (format == "spectronaut") {

      spectraLib[, c("FragmentMz", "FragmentCharge", "RelativeIntensity", "FragmentType", "FragmentNumber")] =
        msLib@PredLib@PrositLib$MsmsData[msmsIdx, c("FragmentMz", "FragmentCharge", "RelativeIntensity", "FragmentType", "FragmentNumber")]
      spectraLib$ModifiedPeptide = str_c("_", str_replace_all(spectraLib$StrippedPeptide, "C", "C[Carbamidomethyl (C)]"), "_")
      spectraLib$FragmentLossType = "noloss"
      del = ","

    } else if (format == "openswath") {

      spectraLib[, c("ProductMz", "ProductCharge", "LibraryIntensity", "FragmentType", "FragmentSeriesNumber")] =
        msLib@PredLib@PrositLib$MsmsData[msmsIdx, c("FragmentMz", "FragmentCharge", "RelativeIntensity", "FragmentType", "FragmentNumber")]
      print("Adding OpenSwath Grouping Columns to TSV Library...")
      spectraLib[, c("DetectingTransition", "QuantifyingTransition")] = 1
      spectraLib[, c("IdentifyingTransition", "Decoy")] = 0
      spectraLib[, c("CollisionEnergy", "PrecursorIonMobility")] = -1
      spectraLib$TransitionId = 1:nrow(spectraLib)
      spectraLib$ModifiedPeptideSequence = str_replace_all(spectraLib$PeptideSequence, "C", "C(UniMod:4)")
      spectraLib$TransitionGroupId = str_c(str_replace_all(spectraLib$PeptideSequence, "C", "C(Carbamidomethyl)"), "_", spectraLib$PrecursorCharge)
      spectraLib$Annotation = str_c(spectraLib$FragmentType, spectraLib$FragmentSeriesNumber, "^", spectraLib$ProductCharge)
      del = "\t"
    }

    if(i == 1) {
      append = F
    } else {
      append = T
    }
    print(str_c("Writing batch ", as.character(i), "/", as.character(nbrOfAppends), " to ", outputLib))
    write_delim(x = spectraLib, file = outputLib, delim = del, append = append) # Is not true for spectronaut...
  }
  outputLib
}


