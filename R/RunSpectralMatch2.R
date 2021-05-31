#' Runs spectral comparisons between experimental spectra in the calibration library object and predicted spectra
#' @param spectra input Spectrum2 objects of Calibration Library object
#' @param precursors filtered precursor data in Calibration library object
#' @param predictionDb absolute path to the prediction database
#' @param threads integer giving the number of threads to use for computation
#' @export run.spectral.match2

run.spectral.match2 <- function(spectra, precursors, predictionDb, threads) {

  add.msms.data <- function(count, precursorData, predSpectra, prositdb, ce) {

    queryStr = str_c("SELECT RelativeIntensity_",ce, " FROM MsmsData WHERE rowid >= ",
                     as.character(precursorData$startIdx[count]), " AND rowid <=",
                     as.character(precursorData$endIdx[count]))
    predSpectra[[count]]@intensity = dbGetQuery(prositdb, queryStr)[,1]
    queryStr = str_c("SELECT FragmentMz FROM MsmsData WHERE rowid >= ",
                     as.character(precursorData$startIdx[count]), " AND rowid <=",
                     as.character(precursorData$endIdx[count]))
    predSpectra[[count]]@mz = dbGetQuery(prositdb, queryStr)[,1]
    predSpectra[[count]]@intensity = predSpectra[[count]]@intensity[order(predSpectra[[count]]@mz)]
    predSpectra[[count]]@mz = predSpectra[[count]]@mz[order(predSpectra[[count]]@mz)]
    predSpectra[[count]]@peaksCount = length(predSpectra[[count]]@intensity)
    predSpectra[[count]]

  }
  compare.spectra <- function(count, exp, pred) {

    print(count)
    dotp = compareSpectra(trimMz(exp[[count]], c(200, Inf)),
                          trimMz(pred[[count]], c(200, Inf)), fun = "dotproduct", 0.02)

  }
  #threads = 12
  #predictionDb = "D:/Data_PROSIT/Libraries/Human_Plasmo/SQLITE/human_plasmo_merged_prositdb.sqlite"
  print("Adding precursor data in the prediction database to Spectrum2 objects...")
  prositdb <- dbConnect(RSQLite::SQLite(), predictionDb)
  precursorData = dbGetQuery(prositdb, "SELECT LabeledPeptide, PrecursorCharge, startIdx, endIdx FROM PrecursorData") # Extract PrecursorData from Prosit DB
  intCols = colnames(dbGetQuery(prositdb, "SELECT * FROM MsmsData LIMIT 1"))
  intCols = intCols[grep("RelativeIntensity", intCols)]
  ceSel = sort(as.numeric(str_extract(intCols, pattern = "\\d+")))
  dbDisconnect(prositdb)
  rownames(precursorData) = str_c(precursorData$LabeledPeptide, "_", precursorData$PrecursorCharge)
  precursorData[,c("LabeledPeptide", "PrecursorCharge")] = NULL
  precursorData$precIdx = seq(1, nrow(precursorData))
  precursorStr = str_c(precursors$PeptideSequence,
                     "_",
                     precursors$PrecursorCharge)
  precursorData = precursorData[precursorStr, ]
  gc()
  if(identical(precursorStr[!is.na(precursorData$precIdx)], rownames(precursorData)[!is.na(precursorData$precIdx)])) {

    print("Precursors match...")
    # spectra = calibrationLib@Spectra[calibrationLib@PrecursorData$FilterLib$msmsIdx]
    spectra = spectra[precursors$msmsIdx]
    spectra = spectra[!is.na(precursorData$precIdx)]
    predSpectra = spectra
    precursorData = precursorData[!is.na(precursorData$precIdx),]
    dotpMat = matrix(NA, nrow = nrow(precursorData), ncol = length(ceSel))
    colnames(dotpMat) = str_c("ce_", ceSel)
    gc()
    for (i in ceSel) {
      tic()
      print(str_c("Comparing experimental spectra and predicted spectra having  CE = ", i))
      prositdb <- dbConnect(RSQLite::SQLite(), predictionDb)
      predSpectra = lapply(seq(1,nrow(precursorData)), add.msms.data, precursorData, predSpectra, prositdb, i)
      dbDisconnect(prositdb)

      cl = makeCluster(threads)
      clusterEvalQ(cl, {
        .libPaths(.libPaths()) # Sets the library path
        library(MSnbase)
      })
      dotpMat[,str_c("ce_",i)] = unlist(parLapply(cl, 1:length(spectra), compare.spectra, spectra, predSpectra))
      stopCluster(cl)
      toc()
    }
  } else {
    stop("Precursors do not match. Terminating code execution...")
  }
  gc()
  ceMat = as.data.frame(do.call('rbind', str_split(rownames(precursorData), "_")))
  colnames(ceMat) = c("sequence", "charge")
  ceMat = ceMat[rep(seq_len(nrow(precursorData)), length(ceSel)),]
  ceMat$ce = rep(ceSel, each = nrow(precursorData))
  ceMat$dotproduct = as.vector(dotpMat)
  ceMat$peptideLength = nchar(ceMat$sequence)
  rm(dotpMat, precursorData, spectra, predSpectra, intCols, ceSel)
  gc()
  print("Adding spectral matching result to Calibration Library object, slot: Comparisons")
  ceMat

}







