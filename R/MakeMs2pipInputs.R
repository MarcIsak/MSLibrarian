#' Creates all the necessary input files for prediction of fragment ion intensities and retention times in MS2PIP (and DeepLC)
#' @param msLib an input MSLibrarian class object
#' @param model character specifying which deeplearning model to use. Available models are c("HCD", "CID", "iTRAQ", "iTRAQphospho", "TMT", "TTOF5600", "HCDch2", "CID")
#' @param fragError Fragment ion mass tolerance in Da (only relevant if an MGF file is added)
#' @param libFormat character or vector with the desired library file output types, can be c("csv", "mgf", "msp", "spectronaut", "bibliospec")
#' @param predictData character indicating if precursors should be predicted from the Sequence object(custom) or the Spectrast object(experimental)
#' @param predictRT Logical indicating whether input files should be generated, for retention time predictions in DeepLC.
#' @param sampleNbr Integer specifying a random sample number for retention time calibration.
#' @param batchSize Integer setting the maximum number of entries of a generated peprec batch file
#' @param outputFolder Output folder where all files will be added
#' @param threads integer specifying the number of threads to use for computation
#' @param outputFolder character specifying the absolute path to the folder where all output files will be added.
#' @export make.ms2pip.inputs

make.ms2pip.inputs <- function(msLib, model, fragError, libFormat, predictData, predictRT, sampleNbr, batchSize, threads, outputFolder) {

  unique.prec <- function(charge, prec) {

    print(charge)
    tmp = prec[as.numeric(prec[,"Charge"]) == charge,]
    tmp[!duplicated(tmp[,"Sequence"]),]

  }

  as.peprec <- function(x) {

    counter <<- counter + 1
    mods = ""
    seq = as.character(x[1])

    if(str_count(seq, "ox") != 0) {

      seq = as.character(x[1])
      i = 1
      oxPos = NA

      while(str_count(seq, "ox") != 0) {

        oxPos[i] = str_locate(seq, "M[(]ox[)]")[[1]]
        seq = str_replace(seq,"[(]ox[)]","")
        i = i + 1

      }
      mods = paste(as.character(oxPos), "|","Oxidation","|",sep = "",collapse ="")

    }
    if(str_count(seq, "C") != 0) {

      carbPos = str_locate_all(seq,"C")[[1]][,"start"]
      carbMods = paste(as.character(carbPos), "|","Carbamidomethyl","|",sep = "",collapse ="")
      mods = paste(mods, carbMods, sep = "", collapse = "")

    }

    if(str_count(seq, "ox") == 0 & str_count(seq, "C") == 0) {
      mods = "-|"

    }

    c(paste("peptide", as.character(counter), sep = ""),
          substr(mods,1,nchar(mods)-1),
          seq,  # Removes the last '|' symbol from the modified string.
          as.character(x[2]))
  }

  if(predictData == "Sequence") {

    inPrec = do.call('rbind',lapply(sort(unique(as.numeric(msLib@Sequences@Precursors[,"Charge"]))),
                                    unique.prec, msLib@Sequences@Precursors))
    preCols = c("Sequence", "Charge")
  } else if(predictData == "Spectrast") {
    inPrec = msLib@SpectrastLib@Unfiltered@PrecursorData$FilterLib
    preCols = c("ModifiedPeptideSequence", "PrecursorCharge")
    print("got into the else if")
  }

  print("Creating PEPREC file")
  counter <<- 0
  peprec = t(apply(inPrec[,preCols], 1, as.peprec))
  colnames(peprec) = c("spec_id","modifications", "peptide", "charge")
  print("Writing the PEPREC file...")

  fileNbr = ceiling(nrow(peprec)/batchSize)
  residual = nrow(peprec) - floor(nrow(peprec)/batchSize)*batchSize
  fileIdx = 0

  for (i in 1:fileNbr) {

    if(i == fileNbr) {
      wFile = peprec[(1:residual + fileIdx),]
    } else {
      wFile = peprec[(1:batchSize + fileIdx),]
      fileIdx = batchSize*i
    }
    write.table(wFile,file=paste(outputFolder, "precursors_",as.character(i),".peprec", sep = ""),
                append=FALSE, quote = FALSE, sep = ",", row.names = F, col.names = T)
  }

  config = c(paste("model=",model,sep=""),
             paste("frag_error=",fragError, sep=""),
             paste("out=", libFormat, sep=""),
             "ptm=Carbamidomethyl,57.02146,opt,C")
  print("Writing the config.txt file...")
  write.table(config,file=paste(outputFolder, "config_ms2pip.txt", sep = ""),append=FALSE, quote = FALSE, sep = "\t", row.names = F, col.names = F)

  if(predictRT == T & predictData == "Sequence") {

    print("Creating input files for MS2PIP RT predictions...")
    inPrec = msLib@SpectrastLib@Unfiltered@PrecursorData$FilterLib
    counter <<- 0
    calib = t(apply(inPrec[,c("ModifiedPeptideSequence", "PrecursorCharge")],1, as.peprec))
    colnames(calib) = c("spec_id","modifications", "peptide", "charge")
    selInPrec = inPrec[!duplicated(calib[,"peptide"]),] # temp added in here....
    calib = calib[!duplicated(calib[,"peptide"]),]
    #selInPrec = inPrec[!duplicated(calib[,"peptide"]),] # temporarily added out...
    t = sample(nrow(calib), sampleNbr, replace = F)
    calibrSet = cbind(calib[t, c("peptide","modifications")],selInPrec[t,"NormalizedRetentionTime"])
    # Something is going wrong here when adding retentiontimes to the calib dataset....

  } else if (predictRT == T & predictData == "Spectrast") {

    selPepRec = peprec[!duplicated(peprec[,"peptide"]),]
    selInPrec = inPrec[!duplicated(peprec[,"peptide"]),]
    t = sample(nrow(selPepRec), sampleNbr, replace = F)
    calibrSet = cbind(selPepRec[t, c("peptide","modifications")],selInPrec[t,"NormalizedRetentionTime"])
  }

    calibrSet[,"modifications"] = str_remove(calibrSet[,"modifications"],"-")
    colnames(calibrSet) = c("seq", "modifications", "tr")
    #predSet = peprec[!duplicated(peprec[,"peptide"]), c("peptide","modifications")]
    predSet = peprec[, c("peptide","modifications")]
    colnames(predSet) = c("seq", "modifications")
    predSet[, "modifications"] = str_remove(predSet[, "modifications"],"-")
    print("Write RT prediction input files...")

    write.table(predSet, file=paste(outputFolder, "precursors","_rt.csv", sep = ""),append=FALSE, quote = FALSE, sep = ",", row.names = F, col.names = T)
    write.table(calibrSet, file=paste(outputFolder, "precursors","_rt_calib.csv", sep = ""),append=FALSE, quote = FALSE, sep = ",", row.names = F, col.names = T)

}

