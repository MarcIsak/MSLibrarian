#' Draws a scatter plot of either dotproducts or pearsons correlation coefficients for spectral matches.
#' @param pepShkr input PeptideShaker object
#' @param mode character specifying if dot products or pearsons correlation coefficient should be plotted.
#' @param matches character specifying which spectral matches to plot against each other(x vs. y)
#' @param xlab character specifying the x-axis label
#' @param ylab character specifying the y-axis label
#' @param main character specifying the title of the plot
#' @export scatter.spectral.match

scatter.spectral.match <- function(pepShkr,mode,matches, xlab, ylab, main) {

  plot(pepShkr@psmsPass[[matches[1]]][,mode][!is.na(pepShkr@psmsPass[[matches[1]]][,mode])],
       pepShkr@psmsPass[[matches[2]]][,mode][!is.na(pepShkr@psmsPass[[matches[2]]][,mode])],
       xlab = xlab, ylab = ylab, main = main,
       xlim = c(0,1),ylim = c(0,1),cex = 1.2,font.lab = 2)
  abline(0,1, col = "red")


}
