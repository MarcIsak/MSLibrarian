#' Performs a filtering of MSLibrarian-predicted libraries based on extracted features
#' @param inputLib input spectral library in Spectronaut .csv format
#' @param outputLib filtered output spectral library in Spectronaut .csv format
#' @param dinoPath path to the Dinosaur application (.jar) for feature extraction.
#' @param mzml path to a mzml file. If a features file is supplied, this parameter is ignored.
#' @param tolMz integer specifying the Precursor M/Z tolerance in PPM (default:10)
#' @param tolRt allowed retention time tolerance in minutes for targets in library (default:5)
#' @param threads number of threads to use for feature matching
#' @export feature.filtering

feature.filtering <- function(inputLib, outputLib, dinoPath, mzml, tolMz, tolRt, threads) {


  feature.pass <- function(mat, features) {

    any(dplyr::between(features$mostAbundantMz[features$charge == mat["PrecursorCharge"]],
                       left = mat["lowMz"],
                       right = mat["highMz"]) &
          dplyr::between(features$rtApex[features$charge == mat["PrecursorCharge"]],
                         left = mat["lowRt"],
                         right = mat["highRt"]))

  }
  pass.msms.idx <- function(mat) {
    seq(mat["startIdx"], mat["endIdx"])

  }
  print("Importing library...")
  specLib = read_csv(inputLib)
  print("Extract precursor information...")
  startIdx = which(specLib$FragmentType == "y" & specLib$FragmentNumber == 1 & specLib$FragmentCharge == 1)
  precursorData = specLib[startIdx,c("PrecursorMz", "PrecursorCharge", "iRT")]
  precursorData$endIdx = c(startIdx[2:length(startIdx)]-1, nrow(specLib))
  precursorData$startIdx = startIdx
  print(str_c("Found ", nrow(precursorData), " precursors..."))
  chargeRange = unique(specLib$PrecursorCharge)
  rm(startIdx, specLib)
  gc()
  print(str_c("Determine experimental precursor m/z ranges for: ", tolMz, " ppm..."))
  precursorData$lowMz =  precursorData$PrecursorMz*(1 - tolMz*10^-6)
  precursorData$highMz = precursorData$PrecursorMz*(1 + tolMz*10^-6)
  print(str_c("Determine experimental rt ranges for set RT tolerance: ", tolRt, " min..."))
  precursorData$lowRt = precursorData$iRT - tolRt
  precursorData$highRt = precursorData$iRT + tolRt
  precursorData$negativeRt = precursorData$iRT <= 0
  passIdx = vector("list", length(mzml))

  for (i in 1:length(mzml)) {
    if(file.exists(str_replace(mzml[i], pattern = ".mzML$", ".features.tsv"))) {
      print(str_c("Found feature file for file: ", mzml[i], "..."))
      print("Will use this feature file for further processing...")
      features = str_replace(mzml[i], pattern = ".mzML$", ".features.tsv")
      print(str_c("Importing feature file from: ", features))
      features = read_tsv(features)[,c("mostAbundantMz", "charge", "rtStart", "rtEnd", "rtApex")]
      gc()
    } else if(file.exists(mzml[i]) & str_detect(mzml[i], ".mzML$") & length(mzml[i]) == 1){

      cat("No feature file supplied.\nRunning Dinosaur on: " , mzml[i]) # Perhaps make this into an exported function that is called
      arguments = c("-jar",
                    dinoPath,
                    str_c("--minCharge=", min(chargeRange)),
                    str_c("--maxCharge=", max(chargeRange)),
                    "--verbose",
                    mzml[i])
      system2(command = "java",
              args = arguments,
              wait = T,
              invisible = F,
              minimized = F)
      Sys.sleep(2)
      print("Importing feature file...")
      features = read_tsv(str_replace(mzml[i], pattern = ".mzML$", ".features.tsv"))[,c("mostAbundantMz", "charge", "rtStart", "rtEnd", "rtApex")]
      gc()
    } else {
      print("Incorrect or no files supplied...terminating code execution")
      stop()
    }

    print("Find library entries that match extracted features...")
    cl = parallel::makeCluster(threads)
    clusterEvalQ(cl,{
      .libPaths("C:/Users/marc/Dropbox (Human Neural Develop)/Marc/MISTR/R/win-library/3.5")
      library(MSLibrarian)
    })
    passIdx[[i]] = which(parallel::parApply(cl, precursorData, 1, feature.pass, features) & !precursorData$negativeRt)
    print(str_c("Features found for: ", length(passIdx[[i]]), "/", nrow(precursorData)," precursors  (", round(length(passIdx[[i]])*100/nrow(precursorData), 2), " %)"))
    stopCluster(cl)
  }
  totPass = unique(unlist(passIdx))
  print(str_c("Unique features found for: ", length(totPass), "/", nrow(precursorData)," precursors  (", round(length(totPass)*100/nrow(precursorData), 2), " %)"))
  print("Filter library")
  idx = unlist(apply(precursorData[totPass, c("startIdx", "endIdx")], 1, pass.msms.idx))
  rm(precursorData)
  gc()
  print(str_c("Writing library to: ", outputLib))
  write_csv(x = read_csv(inputLib)[idx,], file = outputLib)
  print("Done!")

}
