#' Calculates iRT values for each PSM in psmsPass$SpecLib and plots iRT for predicted vs. experimental
#' @param pepShkr input PeptideShaker object
#' @param prositLib input predicted spectral library
#' @export get.psm.irt

get.psm.irt <- function(pepShkr,prositLib) {

  irtcalc <- function(x) {

    out = pepShkr@psmsPass$iRT_Summary$coefficients[1,"Estimate"] + (pepShkr@psmsPass$iRT_Summary$coefficients[2,"Estimate"]*x)

  }
  get.prosit.irt = function(x) {

    out = rtime(prositLib$Spectrum2[[x]])

  }

  pepShkr@psmsPass$specLib$iRT_calc = sapply(as.numeric(pepShkr@psmsPass$specLib$`retention time`)/60,irtcalc)
  plot(sapply(which(!is.na(pepShkr@psmsPass$specLib$prositIndx)),get.prosit.irt),
       pepShkr@psmsPass$specLib$iRT_calc[!is.na(pepShkr@psmsPass$specLib$prositIndx)],
       ylab = "Predicted iRT : Prosit",xlab = "Experimental iRT") # Plots predicted vs. experimental iRT
  abline(0,1,col="red", lwd = 2)
  pepShkr

}
