#' Plots iRT vs. RT with a linear regression. It also returns a summary object of the linear regression
#' @param pepShkr input PeptideShaker object.
#' @export irt.qc.plot

irt.qc.plot <- function(pepShkr) {

  library("stringr", lib.loc="C:/Users/marc/Dropbox (Human Neural Develop)/Marc/MISTR/R/win-library/3.5")

  get.rt = function(z){

    median(as.numeric(pepShkr@psmsPass$PSMs$`retention time`[which(pepShkr@psmsPass$PSMs$peptide_ref == z)]))
  }

  irtSeqs = c("LGGNEQVTR","GAGSSEPVTGLDAK","VEATFGVDESNAK","YILAGVENSK","TPVISGGPYEYR","TPVITGAPYEYR",
              "DGLDAASYYAPVR", "ADVTPADFSEWSK","GTFIIDPGGVIR","GTFIIDPAAVIR","LFLQFGAQGSPFLK")
  irtIndx = c(-24.92,0.0,12.39,19.79,28.71,33.38,42.26,54.62,70.52,87.23,100.0)
  peptIds = str_replace_all(irtSeqs,"L","I") # Makes sequences into corresponding peptideIds (no L, only I)
  rt = sapply(peptIds,get.rt)/60
  plot(rt,irtIndx,pch = 20,xlab = "RT / min", ylab = "iRT")
  lmSmry = summary(lm(irtIndx ~ rt))
  abline(lmSmry$coefficients[1,"Estimate"],lmSmry$coefficients[2,"Estimate"],col = "red")
  text(max(axTicks(1)),min(axTicks(2)),paste("iRT = ",as.character(round(lmSmry$coefficients[1,"Estimate"],2))," + ",
                    as.character(round(lmSmry$coefficients[2,"Estimate"],2)),"*RT",sep=""))
  text(max(axTicks(1)), min(axTicks(2)) + 5,
       paste("P-value = ", as.character(round(lmSmry$coefficients[2,"Pr(>|t|)"],
                                              as.integer(abs(floor(log10(lmSmry$coefficients[2,"Pr(>|t|)"]))))))))
  text(max(axTicks(1)), min(axTicks(2)) + 10,
       paste("R^2 = ", as.character(round(lmSmry$r.squared,2))))
  pepShkr@psmsPass[["iRT_Summary"]] = lmSmry
  pepShkr




}
