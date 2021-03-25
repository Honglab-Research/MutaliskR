# Housekeeping functions.
#
# Author: Andy Jinseok Lee


#' @title Returns the PCAWG SBS reference signature data
#'
#' @description This function returns the PCAWG single base subsitution
#' (SBS) mutational signatures data.
#'
#' @param sequencing.type Either 'WES' for whole-exome sequencing or
#' 'WGS' for whole-genome sequencing.
#'
#' @return A data.frame of the PCAWG SBS mutational signatures.
#'
#' @export
#' @importFrom utils read.csv
GetPcawgSbsSignaturesData <- function(sequencing.type) {
  if (sequencing.type == "WGS") {
    f <- system.file("extdata/01_Reference_Signatures/01_PCAWG_Signatures",
                     "SigProfiler_genome_SBS_signatures_20200128.csv",
                     package = "MutaliskR")
  } else {
    print("Invalid sequencing.type parameter. It must be either 'WES' or 'WGS'")
    return(NULL)
  }
  df <- read.csv(file = f,
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 dec = ".")
  return(df)
}

#' @title Return the PCAWG SBS reference signature names
#'
#' @description This function returns the PCAWG SBS mutational signatures names.
#'
#' @param sequencing.type Either 'WES' for whole-exome sequencing or
#' 'WGS' for whole-genome sequencing.
#'
#' @return Returns a character vector of all PCAWG SBS mutational signature names.
#'
#' @export
#' @importFrom utils read.csv
GetPcawgSbsSignaturesNames <- function(sequencing.type) {
  df <- GetPcawgSbsSignaturesData(sequencing.type = sequencing.type)
  signatures.names <- as.character(colnames(df)[2:ncol(df)])
  return(signatures.names)
}

#' @title Returns the PCAWG Platinum SBS reference signature data
#'
#' @description This function returns the PCAWG Platinum single base subsitution
#' (SBS) mutational signatures data (from Petljak et al., Cell 2019).
#'
#' @return A data.frame of the PCAWG Platinum SBS mutational signatures.
#'
#' @export
#' @importFrom utils read.csv
GetPcawgPlatinumSbsSignaturesData <- function() {
  f <- system.file("extdata/01_Reference_Signatures/01_PCAWG_Signatures",
                   "PCAWG_Platinum_Signatures_from_Petljak_Cell_2019.csv",
                   package = "MutaliskR")
  df <- read.csv(file = f,
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 dec = ".")
  return(df)
}

#' @title Return the PCAWG Platinum SBS reference signature names
#'
#' @description This function returns the PCAWG Platinum SBS mutational signatures names.
#'
#' @return Returns a character vector of all PCAWG Platinum SBS mutational signature names.
#'
#' @export
#' @importFrom utils read.csv
GetPcawgPlatinumSbsSignaturesNames <- function() {
  df <- GetPcawgPlatinumSbsSignaturesData()
  signatures.names <- as.character(colnames(df)[2:ncol(df)])
  return(signatures.names)
}

#' @title Returns the Petljak SBS reference signature data
#'
#' @description This function returns the Petljak single base subsitution
#' (SBS) mutational signatures data.
#'
#' @return A data.frame of the Petljak SBS mutational signatures.
#'
#' @export
#' @importFrom utils read.csv
GetPetljakSbsSignaturesData <- function() {
  f <- system.file("extdata/01_Reference_Signatures/02_Petljak_Signatures",
                   "Petljak_et_al_Cell_2019_SBS_signatures.csv",
                   package = "MutaliskR")
  df <- read.csv(file = f,
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 dec = ".")
  return(df)
}

#' @title Return the Petljak SBS reference signature names
#'
#' @description This function returns the Petljak SBS mutational signatures names.
#'
#' @return Returns a character vector of all Petljak SBS mutational signature names.
#'
#' @export
#' @importFrom utils read.csv
GetPetljakSbsSignaturesNames <- function(sequencing.type) {
  df <- GetPetljakSbsSignaturesData()
  signatures.names <- as.character(colnames(df)[2:ncol(df)])
  return(signatures.names)
}

#' @title Returns the Kucab SBS reference signature data
#'
#' @description This function returns the Kucab single base subsitution
#' (SBS) mutational signatures data.
#'
#' @return A data.frame of the Kucab SBS mutational signatures.
#'
#' @export
#' @importFrom utils read.csv
GetKucabSbsSignaturesData <- function() {
  f <- system.file("extdata/01_Reference_Signatures/03_Kucab_Signatures",
                   "Kucab_et_al_Cell_2019_SBS_signatures.csv",
                   package = "MutaliskR")
  df <- read.csv(file = f,
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 encoding = "UTF-8",
                 dec = ".")
  return(df)
}

#' @title Return the Kucab SBS reference signature names
#'
#' @description This function returns the Kucab SBS mutational signatures names.
#'
#' @export
#' @importFrom utils read.csv
GetKucabSbsSignaturesNames <- function() {
  df <- GetKucabSbsSignaturesData()
  signatures.names <- as.character(colnames(df)[2:ncol(df)])
  return(signatures.names)
}

#' @title Returns the Pleguezuelos-Manzano SBS reference signature data
#'
#' @description This function returns the Pleguezuelos-Manzano single base subsitution
#' (SBS) mutational signatures data.
#'
#' @return A data.frame of the Petljak Pleguezuelos-Manzano mutational signatures.
#'
#' @export
#' @importFrom utils read.csv
GetPleguezuelosManzanoSbsSignaturesData <- function() {
  f <- system.file("extdata/01_Reference_Signatures/04_Pleguezuelos-Manzano_Signatures",
                   "Pleguezuelos-Manzano_et_al_Nature_2020_SBS_signature.csv",
                   package = "MutaliskR")
  df <- read.csv(file = f,
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 dec = ".")
  return(df)
}

#' @title Return the Pleguezuelos-Manzano SBS reference signature names
#'
#' @description This function returns the Pleguezuelos-Manzano
#' SBS mutational signatures names.
#'
#' @export
#' @importFrom utils read.csv
GetPleguezuelosManzanoSbsSignaturesNames <- function() {
  df <- GetPleguezuelosManzanoSbsSignaturesData()
  signatures.names <- as.character(colnames(df)[2:ncol(df)])
  return(signatures.names)
}

#' @title Returns the PCAWG SBS penta-nucleotide reference signature data
#'
#' @description This function returns the PCAWG single base substitution
#' (SBS) pentanucleotide context mutational signatures data.
#'
#' @return A data.frame of the PCAWG SBS penta-nucleotide mutational signatures.
#'
#' @export
#' @importFrom utils read.csv
GetPcawgSbsPentanucleotideSignaturesData <- function() {
  f <- system.file("extdata/01_Reference_Signatures/01_PCAWG_Signatures",
                   "SignatureAnalyzer_COMPOSITE_SBS_Pentanucleotide_signatures_031918.csv",
                   package = "MutaliskR")
  df <- read.csv(file = f,
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 dec = ".")
  return(df)
}

#' @title Return the PCAWG SBS penta-nucleotide reference signature names
#'
#' @description This function returns the PCAWG SBS penta-nucleotide mutational signatures names.
#'
#' @return Returns a character vector of all PCAWG SBS penta-nucleotide mutational signature names.
#'
#' @export
#' @importFrom utils read.csv
GetPcawgSbsPentanucleotideSignaturesNames <- function() {
  df <- GetPcawgSbsPentanucleotideSignaturesData()
  signatures.names <- as.character(colnames(df)[2:ncol(df)])
  return(signatures.names)
}

#' @title Returns the COSMIC v3.1 SBS reference signature data
#'
#' @description This function returns the COSMIC v3.1 single base substitution
#' (SBS) mutational signatures data.
#'
#' @return A data.frame of the COSMIC v3.1 SBS mutational signatures.
#'
#' @export
#' @importFrom utils read.csv
GetCosmicSbsSignaturesData <- function() {
  f <- system.file("extdata/01_Reference_Signatures/05_COSMIC_Signatures",
                   "COSMIC_v3.1_SBS_Signatures.csv",
                   package = "MutaliskR")
  df <- read.csv(file = f,
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 dec = ".")
  return(df)
}

#' @title Return the COSMIC v3.1 SBS reference signature names
#'
#' @description This function returns the COSMIC v3.1 SBS mutational signatures names.
#'
#' @return Returns a character vector of all COSMIC v3.1 SBS mutational signature names.
#'
#' @export
#' @importFrom utils read.csv
GetCosmicSbsSignaturesNames <- function() {
  df <- GetCosmicSbsSignaturesData()
  signatures.names <- as.character(colnames(df)[2:ncol(df)])
  return(signatures.names)
}

#' @title Returns the PCAWG DBS reference signature data
#'
#' @description This function returns the PCAWG doublet base subsitution
#' (DBS) mutational signatures data.
#'
#' @return A data.frame of the PCAWG DBS mutational signatures.
#'
#' @export
#' @importFrom utils read.csv
GetPcawgDbsSignaturesData <- function() {
  f <- system.file("extdata/01_Reference_Signatures/01_PCAWG_Signatures",
                   "SigProfiler_genome_DBS_signatures_20200128.csv",
                   package = "MutaliskR")
  df <- read.csv(file = f,
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 dec = ".")
  return(df)
}

#' @title Return the DBS reference signature names
#'
#' @description This function returns the PCAWG DBS mutational signatures names.
#'
#' @return Returns a character vector of all PCAWG DBS mutational signature names.
#'
#' @export
#' @importFrom utils read.csv
GetPcawgDbsSignaturesNames <- function() {
  df <- GetPcawgDbsSignaturesData()
  signatures.names <- as.character(colnames(df)[2:ncol(df)])
  return(signatures.names)
}

#' @title Returns the COSMIC v3.1 DBS reference signature data
#'
#' @description This function returns the COSMIC v3.1 doublet base substitution
#' (DBS) mutational signatures data.
#'
#' @return A data.frame of the PCAWG DBS mutational signatures.
#'
#' @export
#' @importFrom utils read.csv
GetCosmicDbsSignaturesData <- function() {
  f <- system.file("extdata/01_Reference_Signatures/05_COSMIC_Signatures",
                   "COSMIC_v3.1_DBS_Signatures.csv",
                   package = "MutaliskR")
  df <- read.csv(file = f,
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 dec = ".")
  return(df)
}

#' @title Return the COSMIC v3.1 DBS reference signature names
#'
#' @description This function returns the COSMIC v3.1 DBS mutational signatures names.
#'
#' @return Returns a character vector of all COSMIC v3.1 DBS mutational signature names.
#'
#' @export
#' @importFrom utils read.csv
GetCosmicDbsSignaturesNames <- function() {
  df <- GetCosmicDbsSignaturesData()
  signatures.names <- as.character(colnames(df)[2:ncol(df)])
  return(signatures.names)
}

#' @title Returns the PCAWG indel reference signature data
#'
#' @description This function returns the PCAWG small insertion and deletion
#' (indel) mutational signatures data.
#'
#' @param version Either "SigProfiler" or "SignatureAnalyzer".
#'
#' @return A data.frame of the PCAWG indel mutational signatures.
#'
#' @export
#' @importFrom utils read.csv
GetPcawgIndelSignaturesData <- function(version) {
  if (version == "SigProfiler") {
    f <- system.file("extdata/01_Reference_Signatures/01_PCAWG_Signatures",
                     "SigProfiler_genome_ID_signatures_20200128.csv",
                     package = "MutaliskR")
  }
  if (version == "SignatureAnalyzer") {
    f <- system.file("extdata/01_Reference_Signatures/01_PCAWG_Signatures",
                     "SignatureAnalyzer_ID_signature_012718_20200128.csv",
                     package = "MutaliskR")
  }
  df <- read.csv(file = f,
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 dec = ".")
  return(df)
}

#' @title Returns the PCAWG indel signature names
#'
#' @description This function returns the PCAWG indel mutational signatures names.
#'
#' @param version Either "SigProfiler" or "SignatureAnalyzer".
#'
#' @return Returns a character vector of all PCAWG indel mutational signature names.
#'
#' @export
#' @importFrom utils read.csv
GetPcawgIndelSignaturesNames <- function(version) {
  df <- GetPcawgIndelSignaturesData(version = version)
  signatures.names <- as.character(colnames(df)[2:ncol(df)])
  return(signatures.names)
}

#' @title Returns the Pleguezuelos-Manzano indel reference signature data
#'
#' @description This function returns the Pleguezuelos-Manzano small insertion and deletion
#' (indel) mutational signatures data.
#'
#' @return A data.frame of the Pleguezuelos-Manzano indel mutational signatures.
#'
#' @export
#' @importFrom utils read.csv
GetPleguezuelosManzanoIndelSignaturesData <- function() {
  f <- system.file("extdata/01_Reference_Signatures/04_Pleguezuelos-Manzano_Signatures",
                   "Pleguezuelos-Manzano_et_al_Nature_2020_ID_signature.csv",
                   package = "MutaliskR")
  df <- read.csv(file = f,
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 dec = ".")
  return(df)
}

#' @title Returns the Pleguezuelos-Manzano indel signature names
#'
#' @description This function returns the Pleguezuelos-Manzano indel mutational signatures names.
#'
#' @return Returns a character vector of all Pleguezuelos-Manzano indel mutational signature names.
#'
#' @export
#' @importFrom utils read.csv
GetPleguezuelosManzanoIndelSignaturesNames <- function() {
  df <- GetPleguezuelosManzanoIndelSignaturesData()
  signatures.names <- as.character(colnames(df)[2:ncol(df)])
  return(signatures.names)
}

#' @title Returns the COSMIC v3.1 indel reference signature data
#'
#' @description This function returns the COSMIC v3.1 small insertion and deletion
#' (indel) mutational signatures data.
#'
#' @return A data.frame of the COSMIC v3.1 indel mutational signatures.
#'
#' @export
#' @importFrom utils read.csv
GetCosmicIndelSignaturesData <- function() {
  f <- system.file("extdata/01_Reference_Signatures/05_COSMIC_Signatures",
                   "COSMIC_v3.1_INDEL_Signatures.csv",
                   package = "MutaliskR")
  df <- read.csv(file = f,
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE,
                 dec = ".")
  return(df)
}

#' @title Returns the COSMIC v3.1 indel signature names
#'
#' @description This function returns the COSMIC v3.1 indel mutational signatures names.
#'
#' @return Returns a character vector of all COSMIC v3.1 indel mutational signature names.
#'
#' @export
#' @importFrom utils read.csv
GetCosmicIndelSignaturesNames <- function() {
  df <- GetCosmicIndelSignaturesData()
  signatures.names <- as.character(colnames(df)[2:ncol(df)])
  return(signatures.names)
}

#' @title Returns the SBS signatures colors and etiologies
#'
#' @description This function returns the SBS mutational signatures colors and etiologies.
#'
#' @return A data.frame of the SBS mutational signatures colors and etiologies.
#'
#' @export
#' @importFrom utils read.csv
GetSbsSignaturesColors <- function() {
  f <- system.file("extdata/02_Plotting_Metadata",
                   "01_SBS_Signatures_Colors.csv",
                   package = "MutaliskR")
  df <- read.csv(f,
                 header = TRUE,
                 check.names = FALSE,
                 encoding = "UTF-8",
                 stringsAsFactors = FALSE)
  return(df)
}

#' @title Returns the DBS signatures colors and etiologies
#'
#' @description This function returns the DBS mutational signatures colors and etiologies.
#'
#' @return A data.frame of the DBS mutational signatures colors and etiologies.
#'
#' @export
#' @importFrom utils read.csv
GetDbsSignaturesColors <- function() {
  f <- system.file("extdata/02_Plotting_Metadata",
                   "02_DBS_Signatures_Colors.csv",
                   package = "MutaliskR")
  df <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  return(df)
}

#' @title Returns the indel signatures colors and etiologies
#'
#' @description This function returns the indel mutational signatures colors and etiologies.
#'
#' @return A data.frame of the indel mutational signatures colors and etiologies.
#'
#' @export
#' @importFrom utils read.csv
GetIndelSignaturesColors <- function() {
  f <- system.file("extdata/02_Plotting_Metadata",
                   "03_ID_Signatures_Colors.csv",
                   package = "MutaliskR")
  df <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  return(df)
}

#' @title Returns the hsapiens ensembl data
#'
#' @description This function returns the hsapiens ensembl data.
#'
#' @param version Either "GRCh37" or "GRCh38".
#'
#' @export
#' @importFrom utils read.csv
GetHsapiensEnsemblData <- function(version) {
  if (version == "GRCh38") {
    f <- system.file("extdata/04_Ensembl/H_sapiens",
                     "Ensembl_Genes_100_GRCh38_p13_Mart_Export.txt",
                     package = "MutaliskR")
    df <- read.csv(f, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    names(df)[names(df) == "Gene stable ID"] <- "Ensembl_Gene_ID"
    names(df)[names(df) == "Gene stable ID version"] <- "Ensembl_Gene_ID_Version"
    names(df)[names(df) == "Transcript stable ID"] <- "Ensembl_Transcript_ID"
    names(df)[names(df) == "Transcript stable ID version"] <- "Ensembl_Transcript_ID_Version"
    names(df)[names(df) == "Chromosome/scaffold name"] <- "Chr"
    names(df)[names(df) == "Gene start (bp)"] <- "Gene_Start"
    names(df)[names(df) == "Gene end (bp)"] <- "Gene_End"
    names(df)[names(df) == "Gene name"] <- "Gene_Name"
    names(df)[names(df) == "Transcript start (bp)"] <- "Transcript_Start"
    names(df)[names(df) == "Transcript end (bp)"] <- "Transcript_End"
    names(df)[names(df) == "Transcript length (including UTRs and CDS)"] <- "Transcript_Length"
    names(df)[names(df) == "Transcript type"] <- "Transcript_Type"
    names(df)[names(df) == "NCBI gene (formerly Entrezgene) ID"] <- "Entrez_Gene_ID"
    df$Chr <- paste0("chr", df$Chr)
    df <- df[df$Transcript_Type == "protein_coding",]
    return(df)
  } else if (version == "GRCh37") {
    f <- system.file("extdata/04_Ensembl/H_sapiens",
                     "Ensembl_Genes_100_GRCh37_p13_Mart_Export.txt",
                     package = "MutaliskR")
    df <- read.csv(f, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    names(df)[names(df) == "Gene stable ID"] <- "Ensembl_Gene_ID"
    names(df)[names(df) == "Gene stable ID version"] <- "Ensembl_Gene_ID_Version"
    names(df)[names(df) == "Transcript stable ID"] <- "Ensembl_Transcript_ID"
    names(df)[names(df) == "Transcript stable ID version"] <- "Ensembl_Transcript_ID_Version"
    names(df)[names(df) == "Chromosome/scaffold name"] <- "Chr"
    names(df)[names(df) == "Gene start (bp)"] <- "Gene_Start"
    names(df)[names(df) == "Gene end (bp)"] <- "Gene_End"
    names(df)[names(df) == "Gene name"] <- "Gene_Name"
    names(df)[names(df) == "Transcript start (bp)"] <- "Transcript_Start"
    names(df)[names(df) == "Transcript end (bp)"] <- "Transcript_End"
    names(df)[names(df) == "Transcript length (including UTRs and CDS)"] <- "Transcript_Length"
    names(df)[names(df) == "Transcript type"] <- "Transcript_Type"
    names(df)[names(df) == "EntrezGene ID"] <- "Entrez_Gene_ID"
    df$Chr <- paste0("chr", df$Chr)
    df <- df[df$Transcript_Type == "protein_coding",]
    return(df)
  } else {
    print("Invalid version parameter. It must be either 'GRCh37' or 'GRCh38'")
    return(NULL)
  }
}

#' @title Returns the mmusculus ensembl data
#'
#' @description This function returns the mmusculus ensembl data.
#'
#' @param version Either "NCBIM37" or "GRCm38".
#'
#' @export
#' @importFrom utils read.csv
GetMmusculusEnsemblData <- function(version) {
  if (version == "GRCm38") {
    f <- system.file("extdata/04_Ensembl/M_musculus",
                     "Ensembl_Genes_100_GRCm38_p6_Mart_Export.txt",
                     package = "MutaliskR")
    df <- read.csv(f, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    names(df)[names(df) == "Gene stable ID"] <- "Ensembl_Gene_ID"
    names(df)[names(df) == "Gene stable ID version"] <- "Ensembl_Gene_ID_Version"
    names(df)[names(df) == "Transcript stable ID"] <- "Ensembl_Transcript_ID"
    names(df)[names(df) == "Transcript stable ID version"] <- "Ensembl_Transcript_ID_Version"
    names(df)[names(df) == "Chromosome/scaffold name"] <- "Chr"
    names(df)[names(df) == "Gene start (bp)"] <- "Gene_Start"
    names(df)[names(df) == "Gene end (bp)"] <- "Gene_End"
    names(df)[names(df) == "Gene name"] <- "Gene_Name"
    names(df)[names(df) == "Transcript start (bp)"] <- "Transcript_Start"
    names(df)[names(df) == "Transcript end (bp)"] <- "Transcript_End"
    names(df)[names(df) == "Transcript length (including UTRs and CDS)"] <- "Transcript_Length"
    names(df)[names(df) == "Transcript type"] <- "Transcript_Type"
    names(df)[names(df) == "NCBI gene (formerly Entrezgene) ID"] <- "Entrez_Gene_ID"
    df$Chr <- paste0("chr", df$Chr)
    df <- df[df$Transcript_Type == "protein_coding",]
    return(df)
  } else if (version == "NCBIM37") {
    f <- system.file("extdata/04_Ensembl/M_musculus",
                     "Ensembl_Genes_67_NCBIM37_Mart_Export.txt",
                     package = "MutaliskR")
    df <- read.csv(f, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    names(df)[names(df) == "Ensembl Gene ID"] <- "Ensembl_Gene_ID"
    names(df)[names(df) == "Ensembl Transcript ID"] <- "Ensembl_Transcript_ID"
    names(df)[names(df) == "Chromosome Name"] <- "Chr"
    names(df)[names(df) == "Gene Start (bp)"] <- "Gene_Start"
    names(df)[names(df) == "Gene End (bp)"] <- "Gene_End"
    names(df)[names(df) == "Transcript Start (bp)"] <- "Transcript_Start"
    names(df)[names(df) == "Transcript End (bp)"] <- "Transcript_End"
    names(df)[names(df) == "Transcript Biotype"] <- "Transcript_Type"
    names(df)[names(df) == "Associated Gene Name"] <- "Gene_Name"
    names(df)[names(df) == "EntrezGene ID"] <- "Entrez_Gene_ID"
    df$Chr <- paste0("chr", df$Chr)
    df <- df[df$Transcript_Type == "protein_coding",]
    return(df)
  } else {
    print("Invalid version parameter. It must be either 'GRCh37' or 'GRCh38'")
    return(NULL)
  }
}
