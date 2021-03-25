# VCF parsing functions.
#
# Author: Andy Jinseok Lee


#' @title Prepares a VCF file
#'
#' @description This function prepares a VCF file.
#'
#' @param vcf.file VCF file including path.
#'
#' @return A data.frame with the following items:
#' \item{Chr}{Chromosome name.}
#' \item{Pos}{Genomic position.}
#' \item{Ref}{Genomic position.}
#' \item{Alt}{Genomic position.}
#'
#' @export
#' @import bedr
PrepareVcfFile <- function(vcf.file) {
  vcf <- read.vcf(x = vcf.file, split.info = FALSE, split.samples = FALSE, verbose = FALSE)
  vcf$vcf$Chr <- vcf$vcf$CHROM
  vcf$vcf$Pos <- vcf$vcf$POS
  vcf$vcf$Ref <- vcf$vcf$REF
  vcf$vcf$Alt <- vcf$vcf$ALT

  # 'chr' formatting
  if (str_detect(string = vcf$vcf$Chr[1], pattern = "chr") == FALSE) {
    vcf$vcf$Chr <- paste0("chr", vcf$vcf$Chr)
  }
  vcf$vcf[vcf$vcf$Chr == "chrMT", "Chr"] <- "chrM"

  df <- vcf$vcf[,c("Chr", "Pos", "Ref", "Alt")]
  return(df)
}
