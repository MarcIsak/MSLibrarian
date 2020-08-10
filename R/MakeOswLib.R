#' Creates a generic spectral library from predicted data, which can be used in OpenSWATH (TSV-format)
#' @param msLib input MSLibrarian object
#' @param selLib character setting the spectral library format --> c("PrositLib", "Ms2pipLib")
#' @param selRT character defining which indexed retention time predictions to use c("prosit", "ms2pip")
#' @param idSource character specifying where the library comes from c("msLibrarian", "spectrast")
#' @param libOutput character specifying the output library file name (incl. absolute path)
#' @param batchSize numeric setting the maximum number of entries to process at a time (to keep the memory usage low)
#' @export make.osw.lib

make.osw.lib <- function(msLib, selLib, selRT, idSource,libOutput, batchSize) {

  rep.entries <- function(x, selCol) {

    rep(as.character(x[selCol]),
        length(seq(as.numeric(x["startIdx"]),as.numeric(x["endIdx"]))))

  }
  specLibCols = c("PrecursorMz", "ProductMz", "PrecursorCharge", "ProductCharge", "LibraryIntensity", "NormalizedRetentionTime",
                  "PeptideSequence", "ModifiedPeptideSequence", "PeptideGroupLabel", "LabelType", "CompoundName", "SumFormula",
                  "SMILES", "ProteinId", "UniprotId", "FragmentType", "FragmentSeriesNumber", "Annotation", "CollisionEnergy",
                  "PrecursorIonMobility", "TransitionGroupId", "TransitionId", "Decoy", "DetectingTransition", "IdentifyingTransition",
                  "QuantifyingTransition", "Peptidoforms")
  predLibCols = c("PrecursorMz", "PrecursorCharge", "iRT", "StrippedPeptide", "LabeledPeptide", "ID")
  selSpecLibCols = c("PrecursorMz", "PrecursorCharge", "NormalizedRetentionTime", "PeptideSequence",
                     "ModifiedPeptideSequence", "ProteinId")

  if(selLib == "PrositLib" & selRT == "ms2pip" | selLib == "Ms2pipLib" & selRT == "prosit") {
    if(!all(msLib@PredLib@PrositLib$PrecursorData$StrippedPeptide == msLib@PredLib@Ms2pipLib$PrecursorData$StrippedPeptide)) {
      print("Error: Predicted libraries do not match!")
      # Must use the modified peptides column in the future
      quit()
    } else {
      print("Predicted libraries match!")
    }
  }
  if(selLib == "Ms2pipLib") {
    predLibCols[predLibCols == "LabeledPeptide"] = "ModifiedPeptide"
    print("MS2PIP selected as library")
    if(selRT == "ms2pip") {
      irtCol = msLib@PredLib@Ms2pipLib$PrecursorData[,c("iRT", "startIdx", "endIdx")]
      print("DeepLC iRTs used")
    } else if(selRT == "prosit") {
      irtCol = cbind(msLib@PredLib@PrositLib$PrecursorData$iRT, msLib@PredLib@Ms2pipLib$PrecursorData[,c("startIdx", "endIdx")])
      print("Prosit iRTs used")
    } else {
      print("Incorrect iRT selection!")
      quit()
    }
  } else if(selLib == "PrositLib") {
    print("Prosit selected as library")
    if(selRT == "prosit") {
      irtCol = msLib@PredLib@PrositLib$PrecursorData[,c("iRT", "startIdx", "endIdx")]
      print("Prosit iRTs used")
    } else if(selRT == "ms2pip") {
      irtCol = cbind(msLib@PredLib@Ms2pipLib$PrecursorData$iRT, msLib@PredLib@PrositLib$PrecursorData[,c("startIdx", "endIdx")])
      print("DeepLC iRTs used")
    } else {
      print("Incorrect iRT selection!")
      quit()
    }
  } else {
    print("Selected library does not exist!")
  }
  colnames(irtCol)[1] = "iRT"

  if(nrow(slot(msLib@PredLib,selLib)$PrecursorData) <= batchSize) {
    batchSize = nrow(slot(msLib@PredLib,selLib)$PrecursorData)
    print(paste("Setting batchSize to",nrow(slot(msLib@PredLib,selLib)$PrecursorData)))
    nbrOfAppends = 1
    firstIdx = 1
    lastIdx = batchSize
  } else {
    nbrOfAppends = ceiling(nrow(slot(msLib@PredLib,selLib)$PrecursorData)/batchSize)
    firstIdx = c(1,seq(1, nbrOfAppends-1)*batchSize+1)
    lastIdx = c(seq(1,nbrOfAppends-1)*batchSize, nrow(slot(msLib@PredLib,selLib)$PrecursorData))
  }
  if(idSource == "msLibrarian") {
    idCol = as.data.frame(cbind(msLib@Sequences@Precursors$ID, slot(msLib@PredLib,selLib)$PrecursorData[,c("startIdx", "endIdx")]), stringsAsFactors = F)
  } else if(idSource == "spectrast") {
    idCol = as.data.frame(cbind(msLib@SpectrastLib@Unfiltered@PrecursorData$FilterLib$ProteinId,
                                slot(msLib@PredLib,selLib)$PrecursorData[,c("startIdx", "endIdx")]), stringsAsFactors = F)
  } else {
    print("Incorrect ID Source!")
  }
  colnames(idCol)[1] = "ID"
  file.create(libOutput)
  write.table(t(specLibCols), file = libOutput, sep = "\t", col.names = F, row.names = F, quote = F, append = T)
  for(i in 1:nbrOfAppends) {

    firstMsMsIdx = slot(msLib@PredLib,selLib)$PrecursorData$startIdx[firstIdx[i]]
    lastMsMsIdx = slot(msLib@PredLib,selLib)$PrecursorData$endIdx[lastIdx[i]]
    print(firstMsMsIdx)
    print(lastMsMsIdx)

    test = as.data.frame(matrix(NA, length(seq(firstMsMsIdx, lastMsMsIdx)), length(specLibCols)))
    colnames(test) = specLibCols
    print("Dimensions of the subset:")
    print(as.character(dim(test)))
    for (j in 1:length(predLibCols)) {

      print(paste("Creating the", predLibCols[j], "column for subset",as.character(i),"/",as.character(nbrOfAppends),sep = " "))
      if(predLibCols[j] == "ID") {
        test[,selSpecLibCols[j]] = as.vector(apply(idCol[firstIdx[i]:lastIdx[i],], 1, rep.entries, predLibCols[j]))
      } else if(predLibCols[j] == "iRT") {
        test[,selSpecLibCols[j]] = as.vector(apply(irtCol[firstIdx[i]:lastIdx[i],],1, rep.entries, predLibCols[j]))
        print("Adding iRTs")
      } else {
        test[,selSpecLibCols[j]] = as.vector(apply(slot(msLib@PredLib,selLib)$PrecursorData[firstIdx[i]:lastIdx[i],],1, rep.entries, predLibCols[j]))
      }
    }
    if (selLib == "Ms2pipLib") {
      test[,"ModifiedPeptideSequence"] = str_replace_all(test[,"ModifiedPeptideSequence"], "_", "")  # This formatting should be done when importing the library, fix later....
      test[,"ModifiedPeptideSequence"] = str_replace_all(test[,"ModifiedPeptideSequence"], "C\\[\\+57.0\\]", "C(UniMod:4)")
    } else {
      test[,"ModifiedPeptideSequence"] = str_replace_all(test[,"ModifiedPeptideSequence"], "C", "C(UniMod:4)")
    }

    print("Adding fragment ion columns for the subset")
    test[, c("ProductMz",
             "ProductCharge",
             "LibraryIntensity",
             "FragmentType",
             "FragmentSeriesNumber")] = slot(msLib@PredLib,selLib)$MsmsData[firstMsMsIdx:lastMsMsIdx,c("FragmentMz",
                                                                                                       "FragmentCharge",
                                                                                                       "RelativeIntensity",
                                                                                                       "FragmentType",
                                                                                                       "FragmentNumber")]
    print("Adding OpenSWATH Grouping Columns")
    test[,"TransitionGroupId"] = str_c(str_replace_all(test[,"PeptideSequence"], "C", "C(Carbamidomethyl)"), rep("_", nrow(test)), as.character(test[,"PrecursorCharge"]))
    test[, c("CollisionEnergy","PrecursorIonMobility")] = -1
    test[,c("DetectingTransition", "QuantifyingTransition")] = 1
    test[,c("Decoy", "IdentifyingTransition")] = 0
    test[,"TransitionId"] = seq(firstMsMsIdx, lastMsMsIdx)
    test$NormalizedRetentionTime = as.numeric(test$NormalizedRetentionTime)
    test$LibraryIntensity = as.numeric(test$LibraryIntensity)
    print(paste("Writing subset to: ",libOutput, sep = ""))
    write.table(test, file = libOutput, sep = "\t", col.names = F, row.names = F, quote = F, append = T)
    # Consider using write_tsv from tidyverse later as it should be much faster than write.table...
    rm(test)
    gc()
    print("Done writing subset!")

  }
  rm(idCol, irtCol)
  gc()
  toc()
}


