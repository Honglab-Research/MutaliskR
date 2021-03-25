# MutaliskR

```MutaliskR``` is an R software package that enables general-purpose mutational signature identification. The underlying method of MutaliskR is taken from [Mutalisk (PMID 29790943)](https://academic.oup.com/nar/article/46/W1/W102/5001159), which uses Maximum Likelihood Estimation (MLE) to identify mutational signatures in a tumor sample. The original web-based [Mutalisk](http://mutalisk.org) implementation performs signature analysis of single base substitutions (SBS). This R package extends the Mutalisk (MLE) method to facilitate identification of doublet base substitution (DBS) and small insertion and deletion (INDEL) mutational signatures. The variant/gene level signature analysis implemented in this package also uses the method previously published by [Letouze et al. (PMID 29101368)](https://www.nature.com/articles/s41467-017-01358-x).

## 01. Installation
```
install.packages("devtools")
devtools::install_github("Honglab-Research/MutaliskR", dependencies = TRUE)
```

## 02. Dependencies
- R (>= 3.5.0)
- ggplot2
- ggpubr
- svglite
- dplyr
- plyr
- stringr
- bedr
- GenomicRanges
- DNAcopy
- lsa
- scales
- doParallel
- foreach
- parallel
- BSgenome.Hsapiens.UCSC.hg19 or BSgenome.Hsapiens.UCSC.hg38

## 03. Usage

Examples below use ANNOVAR txt files as input. To use VCF files as input, please refer to the example scripts in the ```examples``` folder.

Identification of ```SBS``` Mutational Signatures
```
library(MutaliskR)
library(BSgenome.Hsapiens.UCSC.hg38)

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
  save.dir = "Output directory path"
)
```

Identification of ```DBS``` Mutational Signatures
```
library(MutaliskR)
library(BSgenome.Hsapiens.UCSC.hg38)

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
  save.dir = "Output directory path"
)
```

Identification of ```INDEL``` Mutational Signatures
```
library(MutaliskR)
library(BSgenome.Hsapiens.UCSC.hg38)

annovar.file <- system.file("extdata/03_Samples", "HCC1187--HCC-1187BL.snv.indel.final.v6.annotated.hg38_multianno.txt", package = "MutaliskR")
input <- PrepareAnnovarFile(annovar.file = annovar.file)
results <- IdentifyIndelSignatures(
  input = input,
  bsg = BSgenome.Hsapiens.UCSC.hg38,
  analyze.variants = TRUE,
  analyze.variants.column.gene = "Gene.refGene",
  analyze.variants.column.group = "Func.refGene",
  sample.id = "HCC1187",
  n.cores = 2,
  save.dir = "Output directory path"
)
```

For output examples, please refer to the zipped files in the ```examples``` folder.

## 04. Contact
For any questions or issues, please post them in the Issues tab of this repository.
