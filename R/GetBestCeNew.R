#' Extract the optimal collision energy based on peptide length, peptide length and charge, or just precursor charge
#' @param ceMat input collision energy matrix as a result from running run.spectral.match
#' @export get.best.ce.new

get.best.ce.new <- function(ceMat) {

  get.dotproduct <- function(x, n, mat) {

    #tmp = median(mat$dotproduct[mat$ce == x & mat$peptideLength == n], na.rm = T)
    tmp = median(mat$dotproduct[mat$ce == x & mat$lenRange == n], na.rm = T)

    if(!is.na(tmp)) {
      tmp
    } else {
      tmp = 0
    }
  }
  set.pept.length = function(n, ceSel, mat) {

    dotp = sapply(ceSel, get.dotproduct, n, mat)
    #print(dotp) # remove after
    if(any(dotp != 0)) { # Changed from all...
      k = c(ceSel[dotp == max(dotp)], max(dotp), n)
      #print(k) # remove after
      if(str_detect(n, "_")) {
        n = as.numeric(unlist(str_split(n, pattern = "_")))
        n = seq(min(n), max(n))
        k = t(replicate(n = length(n), expr = k, simplify = T))
        k[,3] = n
        k
        #print(k)
      } else {
        k
      }
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
  print("First")
  ceMat = ceMat[!is.na(ceMat$dotproduct),]
  ceSel = sort(unique(ceMat$ce))

  lenRange = unique(ceMat$lenRange)
  lenRange = lenRange[order(as.numeric(unlist(lapply(str_split(lenRange, "_"), min))))]
  print("Second")
  bestCe = rbind(as.data.frame(do.call('rbind', lapply(lenRange, set.pept.length, ceSel, ceMat[ceMat$charge == 2,]))),
                 as.data.frame(do.call('rbind', lapply(lenRange, set.pept.length, ceSel, ceMat[ceMat$charge == 3,]))))
  print("Third")
  # No hard coding, change how bestCe is calculated...

  colnames(bestCe) = c("ce", "dotproduct", "peptideLength")
  bestCe = bestCe[!is.na(bestCe$dotproduct),]
  bestCe$charge = rep(c(unique(ceMat$charge)),
                      each = length(unique(bestCe$peptideLength)))
  bestCe = do.call('rbind',lapply(unique(bestCe$charge), add.ce, bestCe))
  print("Fourth")
  # Strange messages...look up...

  bestCe = rbind(bestCe, do.call('rbind', lapply(sort(unique(ceMat$charge)), set.charge, ceMat)))

  bestCe$ce = as.numeric(bestCe$ce)
  bestCe$dotproduct = as.numeric(bestCe$dotproduct)
  bestCe$peptideLength = as.numeric(bestCe$peptideLength)
  bestCe

}
