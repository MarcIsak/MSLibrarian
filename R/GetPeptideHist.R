#' Plots a histogram over peptide masses or precursor m/z values.
#' @param sco an input scell object
#' @param dataSort Defines whether masses in slot Peptides or FilterPeptides should be plotted.
#' @param massMode specifies whether peptide masses or precursor m/z values should be plotted.
#' @param bins number of bins to plot
#' @param xlims defines the limits for the x-axis
#' @export get.peptide.hist

get.peptide.hist <- function(sco,massMode,dataSort,bins,xlims) {


  if(massMode == "mass"){

    pass = unlist(sco@Digestion) # Extract logical vector indicating which peptide mass that passed the filter set in function filter.peptides.
    massVals = unlist(slot(sco,dataSort)$Masses[!is.na(slot(sco,dataSort)$Masses)]) # Removes NA list elements before adding list elements to a vector
    hist(table(massVals[!pass],massVals[pass]),breaks = bins, xlab = "Peptide Mass / (Da)", ylab = "Number of Masses", xlim = c(min(xlims),max(xlims))) # plots  histogram

    # THE TABLE ARGUMENTS ARE IN THE WRONG FORMAT, TRY TO FIX THIS LATER...BUT NOTICE THAT IT SHOULD BE A  HISTOGRAM BUT NOT A BARPLOT...

  } else {

    hist(sapply(sco@Precursors[,"M/Z"],as.numeric),breaks = bins, xlab = "Precursor M/Z", ylab = "Number of Precursors", xlim = c(min(xlims),max(xlims)))

  }

}
