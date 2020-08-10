#' Generates a CSV file with precursor information, to be used for prediction of MS/MS spectra in Prosit.
#' @param msLib input MSLibrarian object
#' @param collisionEngy Collision energy to use for the prediction in PROSIT.
#' @param chargeRange Charge ranges considered for precursor ions.
#' @param outputFile full path to the output CSV file
#' @export make.prosit.csv

make.prosit.csv <- function(msLib, collisionEngy, chargeRange, outputFile) {


  selPrec = msLib@Sequences@Precursors[as.numeric(msLib@Sequences@Precursors[,"Charge"]) >= min(chargeRange) &&
                                         as.numeric(msLib@Sequences@Precursors[,"Charge"]) <= max(chargeRange),]

  write.csv(data.frame(modified_sequence = selPrec[,"Sequence"],
                       collision_energy = rep(collisionEngy,times=nrow(selPrec)),
                       precursor_charge = selPrec[,"Charge"]),
            file = outputFile,
            quote = F,
            sep="",row.names = F,
            col.names = T)

}
