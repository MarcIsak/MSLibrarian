#' Adds DeepLC retention times to the predicted Prosit library
#' @param inputLib path to spectral library that should have its retention times replaced.
#' @param outputLib desired path to spectral library with replaced retention time.
#' @param calibLib path to calibration library as created by the MSLibrarian package (.RData)
#' @param nCal decimal number giving the proportion of calibration peptides in calibrationLib to use
#' @param deepLcFolder path to folder with DeepLC result files
#' @export add.deeplc.rt

add.deeplc.rt <- function(inputLib, outputLib, calibLib, nCal, deepLcFolder) {

  rep.rt <- function(mat) {
    rep(mat[1], mat[2])
  }
  if(dir.exists(deepLcFolder)) {
    print(str_c("Loading DeepLC results in folder: ", deepLcFolder))

    bench = read_csv(file.path(deepLcFolder,"deeplc_bench_deeplc_predictions.csv"), col_types = cols())
    calib = read_csv(file.path(deepLcFolder,"deeplc_calib.csv"), col_types = cols())
    libRt = read_csv(file.path(deepLcFolder,"deeplc_lib_deeplc_predictions.csv"),
                     col_types = cols(seq = col_character(),
                                      modifications = col_character(),
                                      predicted_tr = col_double()))
  } else {
    stop("DeepLC results folder or necessary result files do not exist. Stopping code execution!")
  }

  if(identical(bench$seq, calibLib@PrecursorData$FilterLib$PeptideSequence)) {
    rts = data.frame(benchRt = bench$predicted_tr,
                     calibRt = calibLib@PrecursorData$FilterLib$NormalizedRetentionTime)
    model = lm(str_c("benchRt", "~","calibRt"), data = rts)
    ggplot(data = rts, mapping = aes(x = benchRt, y = calibRt)) +
      geom_abline(intercept = 0, slope = 1,color = "black", size = 1) +
      geom_hex(bins = 100) +
      geom_smooth(method = "lm", color = "red", size = 1) +
      geom_text(x = min(rts$benchRt), y = max(rts$calibRt), label = str_c("y = ", round(coef(model)[2],2),
                                                  "x + ",
                                                  round(coef(model)[1],2),
                                                  ", R^2 = ",
                                                  round(summary(model)$r.squared,2)),
                size = 8, color = "red", hjust = "left") +
      # geom_abline(intercept = 0, slope = 1, color = "red", size = 4) +
      #geom_point() +
      xlab("Predicted RT") + # Should be changed to iRT if that is the input...could be extracted from calibLib I guess...
      ylab("Calibration library RT" ) # Should be changed to iRT if that is the input...could be extracted from calibLib I guess...
    print(str_c("Saving retention time benchmark plot to ", file.path(deepLcFolder,"bench.pdf")))
    ggsave(filename = file.path(deepLcFolder,"bench.pdf"))
    # Should write the ggplot object as a PDF file...

  } else {

    stop("Benchmark peptides do not match corresponding calibration library peptides. Stopping code execution.")

  }
  if(str_detect(inputLib, ".csv$")) {
    print(str_c("Importing spectronaut library : ", inputLib))
    delim = ","
    seqCol = "StrippedPeptide"
    rtCol = "iRT"
    modPeptCol = "ModifiedPeptide"

  } else if(str_detect(inputLib, ".tsv$")) {
    print(str_c("Importing OpenSwath library: ", inputLib))
    delim = "\t"
    seqCol = "PeptideSequence"
    rtCol = "NormalizedRetentionTime"
    modPeptCol = "ModifiedPeptideSequence"
  }

  specLib = as.data.frame(read_delim(inputLib, col_types = cols(), delim = delim))
  startIdx = which(!duplicated(str_c(specLib[,modPeptCol], "_", specLib$PrecursorCharge)))
  if(identical(libRt$seq, specLib[startIdx, seqCol])) {
    endIdx = c(startIdx[2:length(startIdx)], nrow(specLib) + 1)
    print("Replacing library retention times with DeepLC predicted retention times...")
    specLib[,rtCol] = unlist(apply(cbind(libRt$predicted_tr, endIdx-startIdx), 1, rep.rt))
    gc()
    print(str_c("Writing new library to: ", outputLib))
    write_delim(x = specLib,file = outputLib, delim = delim)

  } else {
    stop("Sequences do not match. Code execution aborted!")
  }
}
