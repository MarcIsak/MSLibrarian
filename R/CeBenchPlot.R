#' Makes a plot over the best dot products and optimal collision energies based on the spectral matching
#' @param ceMat input collision energy matrix
#' @param outLib path to processed calibration library
#' @export ce.bench.plot

ce.bench.plot <- function(ceMat, outLib) {

  plot.pept.dist <- function(data) {

    data = data[!is.na(data$dotproduct),]
    data = data[data$ce == min(data$ce),]
    data$charge = as.factor(data$charge)

    q = ggplot(data = data, mapping = aes(x = peptideLength, color = charge, fill = charge)) +
      geom_histogram(bins = length(seq(min(data$peptideLength), max(data$peptideLength))), position = "identity", alpha = 0.5) +
      scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
      xlab("Precursor length (Number of amino acids)") + ylab("Number of precursors") +
      ggtitle("Precursor distribution") +
      annotate(geom="text", x=max(data$peptideLength) - 4, y=500, label=as.character(paste("Total:", as.character(nrow(data)))), color="black", size = 7) +
      theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 15, face = "bold"),
            axis.title.y = element_text(size = 15, face = "bold"),
            plot.title = element_text(size = 17, hjust = 0.5, face = "bold.italic", color = "grey30"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 16)) +
      scale_x_continuous(breaks = seq(min(data$peptideLength),max(data$peptideLength),1))
  }
  ce.dotp.plot <- function(data, y, yBreaks) {

    if(y == "dotproduct") {
      title = "Top median dot product between experimental and predicted spectra"
      yLab = "Dot product"
    } else if (y == "ce") {
      title = "Optimal collision energy per precursor length and charge"
      yLab = "Collision energy"
    } else {
      stop("The y - argument is not valid!")
    }

    data = get.best.ce(ceMat)
    data$charge = as.factor(data$charge)

    p = ggplot(data = data[!grepl("only_", data$charge),], mapping = aes_string(x = "peptideLength", y = y, color = "charge")) +
      theme_light() +
      geom_line(size = 2) + scale_x_continuous(breaks = seq(min(data$peptideLength, na.rm = T),max(data$peptideLength, na.rm = T),1)) +
      geom_point(size = 3, color = "black") +
      ggtitle(title) +
      xlab("Precursor length (Number of amino acids)") +
      ylab(yLab) +
      scale_color_manual(values = c("red", "blue", "orange")) +
      theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 15, face = "bold"),
            axis.title.y = element_text(size = 15, face = "bold"),
            plot.title = element_text(size = 17, hjust = 0.5, face = "bold.italic", color = "grey30"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 16)) +
      scale_y_continuous(breaks = yBreaks)

  }
  print("Creating results plot after comparisons of experimental vs. predicted spectra")
  tot = (ce.dotp.plot(data = ceMat, y = "dotproduct", yBreaks = seq(0,1,0.1))/
           ce.dotp.plot(data = ceMat, y = "ce", yBreaks = seq(min(ceMat$ce, na.rm = T),max(ceMat$ce, na.rm = T),2))/
           plot.pept.dist(data = ceMat))
  if(file.exists(outLib) & any(str_detect(outLib, ".tsv$|.RData$"))) {
    print(str_c("Saving plot to:", str_replace(outLib, ".tsv$|.RData$", ".pdf")))
    ggsave(filename = str_replace(outLib, ".tsv$|.RData$", ".pdf"), plot = tot, device = "pdf", width = 10, height = 10, units = "in")
  } else {
    stop("Missing outLib argument or invalid format of outLib argument")
  }


}
