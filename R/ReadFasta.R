#' Reads a fasta file of AA sequences and stores the information in slot Proteins of a MSLibrarian Sequence object.
#' @param file A complete path to a fasta file of AA sequences
#' @export read.fasta

read.fasta = function(file){

  split.fasta.header <- function(x) {

    unlist(str_split(x, " |\\|"))

  }

  get.annot <- function(header, annot) {

    if(any(str_detect(header, annot))) {
      str_remove(header[grepl(annot, header)], annot)
    } else {
      NA
    }

  }
  get.fasta.elements <- function(y) {

    fastaList = data.frame(Accession = y[2],
                           "Entry_name" = y[3],
                           "Protein" = str_c(y[4:(grep("OS=",y)-1)],collapse = " "),
                           "Gene_name" = get.annot(y, "GN="),
                           "Protein_evidence" = get.annot(y, "PE="),
                           "Sequence_version" = get.annot(y, "SV="))
  }

  print(str_c("Importing FASTA database: ", file, "..."))
  # fasta = cbind(do.call(rbind,lapply(sapply(names(as.character(readAAStringSet(file,format = "fasta"))),
  #                                           split.fasta.header),get.fasta.elements)),
  #               Headers = names(as.character(readAAStringSet(file,format = "fasta"))),
  #               Sequences = as.character(readAAStringSet(file,format = "fasta")),row.names = NULL)
  fastaRaw = data.frame(headers = names(as.character(readAAStringSet(file))),
                        sequences = as.character(readAAStringSet(file)))
  if(any(grepl("^DECOY.*_|^REV.*_", fastaRaw$headers, ignore.case = T))) {
    print(str_c("Removing ", sum(!grepl("^DECOY.*_|^REV.*_", fastaRaw$headers, ignore.case = T)), "decoy sequences from database..."))
    fastaRaw = fastaRaw[!grepl("^DECOY.*_|^REV.*_", fastaRaw$headers, ignore.case = T),]
    gc()
  }
  fasta = cbind(do.call(rbind,lapply(sapply(fastaRaw$headers,
                                            split.fasta.header),get.fasta.elements)),
                Headers = fastaRaw$headers,
                Sequences = fastaRaw$sequences,row.names = NULL)
  print(str_c("Database contains: ", nrow(fasta), " proteins..."))
  print(str_c("Adding protein data to slot: Proteins of the MSLibrarian object..."))


   msLib = new("MSLibrarian",
               Sequences = new("Sequence", Proteins = fasta[order(fasta$Accession),]),
               PredLib = new("PredictedLibrary"),
               CalibLib = new("CalibrationLibrary"))







}
