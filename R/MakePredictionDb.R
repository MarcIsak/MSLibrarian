#' Creates a predicted spectral library database
#' @param prediction_folder output folder where the predicted spectral library files will be processed and stored
#' @param task_id text file with input prosit csv file names and their respective task ids
#' @param sqlite preferred name of the SQLite database *.sqlite (It will be added to a subfolder called SQLite in the set prediction_folder)
#' @param download should files in the taskid.txt file be downloaded
#' @param threads number of parallel downloads to perform
#' @export make.prediction.db

make.prediction.db <- function(prediction_folder, task_id, sqlite, download = T,threads) {
  tic()
  prosit.download <- function(id) {
    url = str_c("https://www.proteomicsdb.org/prosit/api/download.xsjs?datasetId=", id["task_id"])
    destFile = str_c(prediction_folder, id["pred_file"],".zip")
    download.file(url = url,
                  mode = "wb",
                  destfile = destFile,
                  cacheOK = F) # this code works and it renames files accordingly...
  }
  task_id = read_delim(task_id, delim = "\t")
  if(download) {
    options(timeout = 180)
    print(str_c("Downloading predicted spectral libraries to: ", prediction_folder))
    apply(task_id, 1, prosit.download)
  }
  # cl = parallel::makeCluster(threads)
  # clusterEvalQ(cl, {
  #   .libPaths(.libPaths())
  #   library(utils)
  #   library(stringr)
  # })
  #doParallel::registerDoParallel(cl)
  #parApply(cl, task_id, 1,prosit.download)

  # foreach(i = 1:nrow(task_id), .packages = "utils") %dopar% {
  #   url = str_c("https://www.proteomicsdb.org/prosit/api/download.xsjs?datasetId=", task_id$task_id[i])
  #   destFile = str_c(prediction_folder, task_id$pred_file[i],".zip")
  #   download.file(url = url,
  #                 mode = "wb",
  #                 destfile = destFile, cacheOK = F) # this code works and it renames files accordingly...
  # }
  # print(getOption('timeout'))
  # stopCluster(cl)
  # gc(full = T)
  # for (i in 1:nrow(task_id)) {
  #   url = str_c("https://www.proteomicsdb.org/prosit/api/download.xsjs?datasetId=", task_id$task_id[i])
  #   destFile = str_c(prediction_folder, task_id$pred_file[i],".zip")
  #   download.file(url = url,
  #                 mode = "wb",
  #                 method = "wininet",
  #                 destfile = destFile) # this code works and it renames files accordingly...
  # }
  ceFiles = list.files(prediction_folder, pattern = "csv.zip")
  ceSel = as.numeric(str_extract(ceFiles, pattern = "\\d+"))
  ceFiles = ceFiles[order(ceSel, decreasing = F)]
  ceSel = ceSel[order(ceSel)]
  toc()
  tic()
  print(str_c("Processing of spectral libraries..."))
  for (i in ceFiles) {

    tic()
    filename = paste(prediction_folder, i, sep = "")
    print(paste("Library to be processed:", filename))
    process.prosit.lib(zipfile = filename, # Could vroom speed up this process perhaps?
                       noZeroIntensity = F,
                       transitionFilter = F,
                       transitionThreshold = 1,
                       output = str_replace(filename, pattern = ".zip", replacement = ".RData"))
    gc(full = T)
    toc()
  }
  toc()
  print(str_c("Merging of spectral libraries into a single database..."))
  ceFiles = str_replace(ceFiles, pattern = ".zip", replacement = ".RData")
  tic()
  for(i in ceFiles) {

    filename = paste(prediction_folder, i, sep = "")
    print(paste("Loading: ", filename, sep = ""))
    load(filename)
    if(i == min(ceFiles)) {
      print("First library loaded!")
      newLib = outLib
      colnames(newLib$MsmsData)[grepl("RelativeIntensity", colnames(outLib$MsmsData))] = paste("RelativeIntensity_",
                                                                                               as.numeric(str_extract(i, pattern = "\\d+")),
                                                                                               sep = "") # Change this...
    } else if (all(newLib$PrecursorData[,c("StrippedPeptide", "startIdx", "endIdx")] ==
                   outLib$PrecursorData[,c("StrippedPeptide", "startIdx", "endIdx")])) {

      newLib$MsmsData[ ,paste("RelativeIntensity_",
                              as.numeric(str_extract(i, pattern = "\\d+")), sep = "")] = outLib$MsmsData$RelativeIntensity # Change this... creates: RelativeIntensity+filename colnames...
      gc(full = T)
    } else {
      print("Sub-libraries do not match! Stopping execution")
      stop()
    }
    rm(outLib)
    gc()
  }
  toc()

  print(str_c("Removal of relative intensities equal to zero..."))
  cols = colnames(newLib$MsmsData)[grep("RelativeIntensity", colnames(newLib$MsmsData))]

  for (i in cols) {

    tmp = which(newLib$MsmsData[,i] == 0)

    if(i == min(cols)) {
      idxRm = tmp
    } else {
      idxRm = intersect(idxRm, tmp)
    }
    print(paste("Done with:", i))
  }

  idxRm = setdiff(idxRm, which(newLib$MsmsData$FragmentType == "y" & newLib$MsmsData$FragmentNumber == 1 & newLib$MsmsData$FragmentCharge == 1))
  # what exists in idxRm that does not exist in the other argument
  newLib$MsmsData = newLib$MsmsData[-idxRm,]
  newLib$PrecursorData$startIdx = which(newLib$MsmsData$FragmentType == "y" & newLib$MsmsData$FragmentNumber == 1 & newLib$MsmsData$FragmentCharge == 1)
  newLib$PrecursorData$endIdx = c(newLib$PrecursorData$startIdx[2:length(newLib$PrecursorData$startIdx)]-1, nrow(newLib$MsmsData))
  gc(full = T)

  print("Saving database as an RData object...")
  save(newLib, file = str_c(prediction_folder,"newLib.RData"), # Temporary to save this object....
       version = NULL,
       ascii = FALSE,
       compress = T,
       safe = T)
  tic()
  sqliteFolder = str_c(prediction_folder, "SQLITE")
  sqlite = str_c(sqliteFolder,"/",sqlite)
  system(str_c("mkdir ", sqliteFolder))
  print(str_c("Exporting the database as an SQLite DB to: ", sqlite))
  prositdb <- dbConnect(RSQLite::SQLite(), sqlite)
  RSQLite::dbWriteTable(conn = prositdb, name = "PrecursorData", value = newLib$PrecursorData)
  RSQLite::dbWriteTable(conn = prositdb, name = "MsmsData", value = newLib$MsmsData)
  dbDisconnect(prositdb)
  toc()

}
