# Examples on doublet base substitution (DBS) mutational signature identification.
#
# Author: Andy Jinseok Lee


library(MutaliskR)
library(BSgenome.Hsapiens.UCSC.hg38)

# Example 1: ANNOVAR (data.frame).
# Variant signatures are analyzed.
annovar.file <- system.file("extdata/03_Samples", "HCC1187--HCC-1187BL.snv.indel.final.v6.annotated.hg38_multianno.txt", package = "MutaliskR")
input <- PrepareAnnovarFile(annovar.file = annovar.file)
results <- IdentifyDbsSignatures(
  input = input,
  bsg = BSgenome.Hsapiens.UCSC.hg38,
  analyze.variants = TRUE,
  analyze.variants.column.gene = "Gene.refGene",
  analyze.variants.column.group = "Func.refGene",
  sample.id = "HCC1187",
  n.cores = 2,
  save.dir = "<PATH>/MutaliskR/examples/02_DBS_Signatures/HCC1187_DBS_Signatures_1/"
)

# Example 2: VCF file.
# Variant signatures are not analyzed.
input <- system.file("extdata/03_Samples", "HCC1187--HCC-1187BL.snv.indel.final.v6.annotated.vcf", package = "MutaliskR")
results <- IdentifyDbsSignatures(
  input = input,
  bsg = BSgenome.Hsapiens.UCSC.hg38,
  sample.id = "HCC1187",
  n.cores = 2,
  save.dir = "<PATH>/MutaliskR/examples/02_DBS_Signatures/HCC1187_DBS_Signatures_2/"
)
