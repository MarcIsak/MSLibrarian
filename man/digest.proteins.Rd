% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DigestProteins.R
\name{digest.proteins}
\alias{digest.proteins}
\title{Returns peptide sequences resulting from protein digestion based on cleavage rules, max.missed cleavages,}
\usage{
digest.proteins(msLib, rowStr, enzyme, maxMissed, carbamidomethyl, threads)
}
\arguments{
\item{msLib}{a MSLibrarian object}

\item{rowStr}{a string defining the row names.}

\item{enzyme}{the enzyme to be used for digestion. The enzyme cleaves proteins according to the specified cleavage rules.}

\item{maxMissed}{maximum missed clevages allowed for a peptide sequence}

\item{carbamidomethyl}{Determines if every Cysteine residue is carbamidomethylated(TRUE) or not (FALSE). If TRUE, every cysteine adds a mass of 57.021464.}

\item{threads}{number of threads to use for the computation}
}
\description{
Returns peptide sequences resulting from protein digestion based on cleavage rules, max.missed cleavages,
}
