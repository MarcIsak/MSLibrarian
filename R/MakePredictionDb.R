#' Creates a predicted spectral library database
#' @param prediction_folder output folder where the predicted spectral library files will be processed and stored
#' @param task_id text file with input prosit csv file names and their respective task ids
#' @param sqlite preferred name of the SQLite database *.sqlite (It will be added to a subfolder called SQLite in the set prediction_folder)
#' @param download should files in the taskid.txt file be downloaded (default = T)
#' @param zipped logical indicating whether downloaded files are zipped (default = T)
#' @export make.prediction.db

make.prediction.db <- function(prediction_folder, task_id, sqlite, download = T, zipped = T) {
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
  if(zipped) {
    files = list.files(prediction_folder, pattern = "csv.zip$", full.names = T)
    ceSel = as.numeric(str_extract(str_extract(basename(files), "_ce.*"), pattern = "\\d+"))
    files = files[order(ceSel, decreasing = F)]
  } else {
    files = list.files(prediction_folder, pattern = ".csv$", full.names = T, include.dirs = T)
    if(!all(files == files[dir.exists(files)])) {
      stop("Not all folders are unzipped.")
    } else if(length(files) == 0) {
      stop("No unzipped folders were found with the ending .csv...")
    } else {
      ceSel = as.numeric(str_extract(str_extract(basename(files), "_ce.*"), pattern = "\\d+"))
      files = files[order(ceSel, decreasing = F)]
      files = file.path(files, "myPrositLib.csv")
    }
  }
  for (i in files) {
    tic()
    print(str_c("Reading file: ", i))
    myPrositLib = read_csv(i, col_types = cols_only(LabeledPeptide = col_character(),
                                                       PrecursorCharge = col_integer(),
                                                       FragmentType = col_character(),
                                                       FragmentNumber = col_integer(),
                                                       FragmentCharge = col_integer()))
    gc()
    if(i == min(files)) {
      precursorComp = str_c(myPrositLib$LabeledPeptide, "_", myPrositLib$PrecursorCharge)
      transitionComp = str_c(myPrositLib$FragmentType, myPrositLib$FragmentNumber, "_", myPrositLib$FragmentCharge)
      firstFile = i
    } else {
      precursors = str_c(myPrositLib$LabeledPeptide, "_", myPrositLib$PrecursorCharge)
      transitions = str_c(myPrositLib$FragmentType, myPrositLib$FragmentNumber, "_", myPrositLib$FragmentCharge)
      if(identical(precursorComp, precursors) & identical(transitionComp, transitions)) {
        print(str_c("Precursors and transitions are identical for ", firstFile, " and ", i))
        rm(precursors, transitions)
        gc(full = T)
      } else {
        rm(myPrositLib, precursors, transitions)
        gc(full = T)
        stop("Precursors and transitions are not identical!")
      }
    }

    rm(myPrositLib)
    gc(full = T)
    toc()
  }
  tic()
  print("Precursor information and transition information agrees for all files...")
  print("Load predicted libraries:")
  for(i in 1:length(files)) { #
    #file = list.files(files[i], pattern = "^myPrositLib.csv$", full.names = T)
    print(str_c("Reading file: ", files[i]))
    if(i == 1) {
      myPrositLib = as.data.frame(read_csv(file = files[i], col_types = cols_only(LabeledPeptide = col_character(),
                                                                                  StrippedPeptide = col_character(),
                                                                                  PrecursorCharge = col_integer(),
                                                                                  PrecursorMz = col_double(),
                                                                                  iRT = col_double(),
                                                                                  FragmentType = col_character(),
                                                                                  FragmentNumber = col_integer(),
                                                                                  FragmentCharge = col_integer(),
                                                                                  FragmentMz = col_double(),
                                                                                  RelativeIntensity = col_double())))
      gc()
    } else {
      print(str_c("Adding predicted fragment ion intensities for Prosit CE: ", ceSel[i]))
      myPrositLib$RelativeIntensity = as.data.frame(read_csv(file = files[i], col_types = cols_only(RelativeIntensity = col_double())))[,1]
      gc()
    }
    colnames(myPrositLib)[grep("^RelativeIntensity$", colnames(myPrositLib))] = str_c("RelativeIntensity_", ceSel[i])
  }
  toc()
  print("Get indices for y1_1 ions...")
  yIonSkip = which(myPrositLib$FragmentType == "y" & myPrositLib$FragmentNumber == 1 & myPrositLib$FragmentCharge == 1)
  for (i in ceSel) {
    tmp = setdiff(which(myPrositLib[,str_c("RelativeIntensity_", i)] == 0), yIonSkip)

    if(i == min(ceSel)) {
      idxRm = tmp
    } else {
      idxRm = intersect(idxRm, tmp)
    }
    print(str_c("CE=", i ,"-> Found ", length(idxRm)," / ", nrow(myPrositLib)," transitions with intensities = 0"))
  }
  print(str_c("Removing ", length(idxRm), " common transitions with intensities = 0 for all collision energies..."))
  myPrositLib = myPrositLib[-idxRm,]
  gc()
  precursorCols = c("LabeledPeptide", "StrippedPeptide", "PrecursorMz", "PrecursorCharge", "iRT")
  msmsCols = setdiff(colnames(myPrositLib), precursorCols)
  startIdx = which(!duplicated(str_c(myPrositLib$LabeledPeptide, "_", myPrositLib$PrecursorCharge)))
  precursorData = myPrositLib[startIdx,precursorCols]
  precursorData$startIdx = startIdx
  precursorData$endIdx = c(startIdx[2:length(startIdx)] - 1, nrow(myPrositLib))
  gc(full = T)
  tic()
  sqliteFolder = file.path(prediction_folder, "SQLITE")
  sqlite = file.path(sqliteFolder, sqlite)
  print(str_c("Creating SQLite folder: ",sqliteFolder))
  dir.create(sqliteFolder)
  print(str_c("Exporting the database as an SQLite DB to: ", sqlite))
  prositdb <- dbConnect(RSQLite::SQLite(), sqlite)
  RSQLite::dbWriteTable(conn = prositdb, name = "PrecursorData", value = precursorData)
  RSQLite::dbWriteTable(conn = prositdb, name = "MsmsData", value = myPrositLib[,msmsCols])
  dbDisconnect(prositdb)
  rm(myPrositLib, precursorData, startIdx, precursorCols, idxRm, msmsCols)
  gc()
  toc()
}
