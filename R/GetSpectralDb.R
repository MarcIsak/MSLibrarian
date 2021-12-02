#' Downloads a Spectral Warehouse SQLite database with predictions of Prosit fragment ion intensities and retention times
#' @param species species for which a SQLite database should be downloaded. Allowed characters c("H.sapiens", "M.musculus", "S.cerevisiae", "E.coli", "C.elegans", "D.melanogaster")
#' @param folder path to download folder. Will default to users Downloads if not specified.
#' @export get.spectral.db

get.spectral.db <- function(species = NULL, folder = NULL) {

  defaultFolder = shell("echo %userprofile%\\downloads", intern = T)

  speciesData = data.frame(species = c("C.elegans",
                                       "H.sapiens",
                                       "D.melanogaster",
                                       "E.coli",
                                       "M.musculus",
                                       "S.cerevisiae"),
                           fileName = c("caenorhabditis_elegans_prosit_intensity_hcd_2020_irt_2019.sqlite",
                                        "human_prosit_intensity_hcd_2020_prosit_irt_2019.sqlite",
                                        "drosophila_melanogaster_prosit_hcd_intensity_2020_irt_2019.sqlite",
                                        "escherichia_coli_k12_prosit_intensity_hcd_prosit_irt.sqlite",
                                        "mouse_prosit_hcd_intensity_2020_irt_2019.sqlite",
                                        "yeast_prosit_hcd_intensity_2020_irt_2019.sqlite"))

  if(is.null(species)) {
    stop("Arg - species: No species provided")
  }
  if(is.null(folder)) {
    message("No download folder provided")
    folder = defaultFolder
  }
  if (dir.exists(folder)) {
    message(cat(str_c("Spectral Warehouse DB database will be downloaded to: " , folder)))
  } else {
    stop("Arg - folder. Provided download folder does not exist.")
  }
  if(!any(species == speciesData$species)) {
    stop("Arg - species: incorrect species specified.")
  } else {
    file = speciesData$fileName[speciesData$species == species]
    message(cat(str_c("Downloading Spectral Warehouse DB for ", species)))
    record = httr::GET(url = "https://zenodo.org/api/records/5749924")
    record = fromJSON(rawToChar(record$content))
    idx = which(record$files$filename == file)
    downloadUrl = record$files$links$download[idx]
    destFile = file.path(folder, file)
    download.file(url = downloadUrl, destfile = destFile)
    message(cat(str_c("Spectral Warehouse DB: ", file, " downloaded to ", folder)))
  }

}
