#' Plots histograms of dot products or pearsons correlation coefficientts from spectral matches (experimental vs. predicted)
#' @param pepShkr input PepShkr object
#' @param mode character specifying whether to plot "dotproduct" or "pearsonsCor"
#' @param matches character vector specifying the spectral matches to include in the histograms.
#' @param titles character vector specifying the titles for each plotted histogram
#' @param xlab character specifying the x-axis label.
#' @param ylim numeric vector specifying the range of the y-axis.
#' @param col character vector specifying the color of the bars
#' @param binSize numeric value giving the number of bins to plot
#' @export hist.spectral.match

hist.spectral.match <- function(pepShkr,mode, matches,titles, xlab, ylim, col, binSize) {

  par(mfrow = c(length(matches),1)) # Should create a custom figure to minimize margins

  for (i in 1:length(matches)) {

    hist(pepShkr@psmsPass[[matches[i]]][,mode],main = titles[i],xlab = xlab,breaks = binSize, col = col[i] ,cex = 1.2,xlim = c(0,1),ylim=ylim)

  }

}
