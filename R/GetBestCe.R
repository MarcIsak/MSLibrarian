#' Extract the optimal collision energy based on peptide length, peptide length and charge, or just precursor charge
#' @param ceMat input collision energy matrix as a result from running run.spectral.match
#' @export get.best.ce

get.best.ce <- function(ceMat) {

  get.dotproduct <- function(x, n, mat) {

    tmp = median(mat$dotproduct[mat$ce == x & mat$peptideLength == n], na.rm = T)
    if(!is.na(tmp)) {
      tmp
    } else {
      tmp = 0
    }
  }
  set.pept.length = function(n, ceSel, mat) {

    dotp = sapply(ceSel, get.dotproduct, n, mat)
    if(any(dotp != 0)) { # Changed from all...
      c(ceSel[dotp == max(dotp)], max(dotp), n)
    } else {
      c(NA, NA, n)
    }
  }
  add.ce <- function(charge, bestCe) {

    ce = bestCe[bestCe$charge == charge,]
    ce$ce[is.na(ce$ce)] = as.numeric(names(sort(table(ce$ce), decreasing = T))[1])
    ce

  }
  get.charge.ce <- function(ce, charge, mat) {

    c(ce, median(mat$dotproduct[mat$ce == ce & mat$charge == charge], na.rm = T), NA, paste("only_",charge, sep = ""))

  }
  set.charge <- function(charge, mat) {

    tmp = t(sapply(sort(unique(mat$ce)), get.charge.ce, charge, mat))
    colnames(tmp) = c("ce", "dotproduct", "peptideLength","charge")
    tmp[tmp[,"dotproduct"] == max(tmp[,"dotproduct"]),]

  }
  #ceMat = calibLib@Comparisons

  ceMat = ceMat[!is.na(ceMat$dotproduct),]
  ceSel = sort(unique(ceMat$ce))

  bestCe = rbind(as.data.frame(do.call('rbind', lapply(sort(unique(ceMat$peptideLength)), set.pept.length, ceSel, ceMat))),
                 as.data.frame(do.call('rbind', lapply(sort(unique(ceMat$peptideLength)), set.pept.length, ceSel, ceMat[ceMat$charge == 2,]))),
                 as.data.frame(do.call('rbind', lapply(sort(unique(ceMat$peptideLength)), set.pept.length, ceSel, ceMat[ceMat$charge == 3,]))))

  # No hard coding, change how bestCe is calculated...

  colnames(bestCe) = c("ce", "dotproduct", "peptideLength")
  bestCe$charge = rep(c("length_only", unique(ceMat$charge)),
                      each = min(which(duplicated(bestCe$peptideLength)))-1)
  bestCe = do.call('rbind',lapply(unique(bestCe$charge), add.ce, bestCe))

  # Strange messages...look up...

  bestCe = rbind(bestCe, do.call('rbind', lapply(sort(unique(ceMat$charge)), set.charge, ceMat)))

  bestCe$ce = as.numeric(bestCe$ce)
  bestCe$dotproduct = as.numeric(bestCe$dotproduct)
  bestCe$peptideLength = as.numeric(bestCe$peptideLength)
  bestCe




}
