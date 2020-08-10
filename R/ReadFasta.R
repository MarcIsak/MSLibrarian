#' Reads a fasta file of AA sequences and stores the information in slot Proteins of a MSLibrarian Sequence object.
#' @param file A complete path to a fasta file of AA sequences
#' @export read.fasta

read.fasta = function(file){

  #library("Biostrings", lib.loc="C:/Users/marc/Dropbox (Human Neural Develop)/Marc/MISTR/R/win-library/3.5") # Make as a dependency later

  split.fasta.header <- function(x) {

    x = gsub("|"," ",as.character(x),fixed = T)
    x = unlist(strsplit(gsub(" ",".",x,fixed=T),".",fixed=T))
    #print(x)

  }
  get.fasta.elements <- function(y) {


    if(any(grepl("GN=",y)) == TRUE) {

      gn = gsub("GN=","",y[grep("GN=",y)],fixed=T)

    } else {
      gn = "NA"
    }
    #print(y)


    fastaList = data.frame(Accession = y[2],
                           "Entry_name" = y[3],
                           "Protein" = paste(y[4:(grep("OS=",y)-1)],collapse = " "),
                           "Gene_name" = gn,
                           "Protein_evidence" = gsub("PE=","",y[grep("PE=",y)],fixed=T),
                           "Sequence_version" = gsub("SV=","",y[grep("SV=",y)],fixed=T))



  }

  fasta = cbind(do.call(rbind,lapply(sapply(names(as.character(readAAStringSet(file,format = "fasta"))),
                                            split.fasta.header),get.fasta.elements)),
                Headers = names(as.character(readAAStringSet(file,format = "fasta"))),
                Sequences = as.character(readAAStringSet(file,format = "fasta")),row.names = NULL)



   msLib = new("MSLibrarian",
               Sequences = new("Sequence", Proteins = fasta),
               PredLib = new("PredictedLibrary"),
               SpectrastLib = new("SpectrastLibrary"),
               Comparisons = list(Consensus = list(),Unfiltered = list()))







}
