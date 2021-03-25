# Examples on single base substitution (SBS) mutational signature identification.
#
# Author: Andy Jinseok Lee


library(MutaliskR)
library(BSgenome.Hsapiens.UCSC.hg38)

# Example 1: ANNOVAR (data.frame).
# Variant signatures are analyzed.
annovar.file <- system.file("extdata/03_Samples", "HCC1187--HCC-1187BL.snv.indel.final.v6.annotated.hg38_multianno.txt", package = "MutaliskR")
input <- PrepareAnnovarFile(annovar.file = annovar.file)
results <- IdentifySbsSignatures(
  input = input,
  bsg = BSgenome.Hsapiens.UCSC.hg38,
  target.signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS6",
                        "SBS8", "SBS10a", "SBS10b", "SBS13", "SBS17a", "SBS17b",
                        "SBS18", "SBS20", "SBS26", "SBS30"),
  analyze.variants = TRUE,
  analyze.variants.column.gene = "Gene.refGene",
  analyze.variants.column.group = "Func.refGene",
  sample.id = "HCC1187",
  n.cores = 2,
  save.dir = "<PATH>/MutaliskR/examples/01_SBS_Signatures/HCC1187_SBS_Signatures_1/"
)

# Example 2: VCF file.
# Variant signatures are not analyzed.
input <- system.file("extdata/03_Samples", "HCC1187--HCC-1187BL.snv.indel.final.v6.annotated.vcf", package = "MutaliskR")
results <- IdentifySbsSignatures(
  input = input,
  bsg = BSgenome.Hsapiens.UCSC.hg38,
  target.signatures = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS6",
                        "SBS8", "SBS10a", "SBS10b", "SBS13", "SBS17a", "SBS17b",
                        "SBS18", "SBS20", "SBS26", "SBS30"),
  sample.id = "HCC1187",
  n.cores = 2,
  save.dir = "<PATH>/MutaliskR/examples/01_SBS_Signatures/HCC1187_SBS_Signatures_2/"
)

# Example 3: Penta-nucleotide context (1536 mutation types).
# Variant signatures are analyzed.
annovar.file <- system.file("extdata/03_Samples", "HCC1187--HCC-1187BL.snv.indel.final.v6.annotated.hg38_multianno.txt", package = "MutaliskR")
input <- PrepareAnnovarFile(annovar.file = annovar.file)
results <- IdentifySbsSignatures(
  input = input,
  bsg = BSgenome.Hsapiens.UCSC.hg38,
  reference = GetPcawgSbsPentanucleotideSignaturesData(),
  target.signatures = GetPcawgSbsPentanucleotideSignaturesNames(),
  context.length = 5,
  analyze.variants = TRUE,
  analyze.variants.column.gene = "Gene.refGene",
  analyze.variants.column.group = "Func.refGene",
  sample.id = "HCC1187",
  n.cores = 2,
  save.dir = "<PATH>/MutaliskR/examples/01_SBS_Signatures/HCC1187_SBS_Signatures_3/"
)
