# Utility functions.
#
# Author: Andy Jinseok Lee


#' @title Returns complement base
#'
#' @description This function returns the complement base of a nucleotide
#'
#' @param nucleotide Nucleotide base.
#'
#' @return "A", "C", "G", or "T"
#'
#' @export
GetComplementBase <- function(nucleotide) {
  if (nucleotide == "T") {
    return("A")
  }
  if (nucleotide == "C") {
    return("G")
  }
  if (nucleotide == "G") {
    return("C")
  }
  if (nucleotide == "A") {
    return("T")
  }
}


#' @title Returns complement bases
#'
#' @description This function returns the complement bases of a nucleotide sequence
#'
#' @param nucleotides Nucleotide bases (e.g. "AG").
#'
#' @return e.g. "TC"
#'
#' @export
GetComplementBases <- function(nucleotides) {
  nucleotides <- strsplit(nucleotides, "")[[1]]
  complements <- rep("N", length(nucleotides))
  for (i in range(1, length(nucleotides))) {
    complements[i] <- GetComplementBase(nucleotides[i])
  }
  return(paste0(complements, collapse = ""))
}

