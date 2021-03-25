# ANNOVAR parsing function.
#
# Author: Andy Jinseok Lee


#' @title Fetches the closest Gene.refGene from an ANNOVAR file
#'
#' @description This function prepares an ANNOVAR file by fetching the closest Gene.refGene from an ANNOVAR file.
#' If the genomic distance information cannot be parsed, the original 'Gene.refGene' value is returned.
#'
#' @param annovar.file ANNOVAR file including path.
#'
#' @return A data.frame with the following items:
#' \item{Chr}{Chromosome name.}
#' \item{Pos}{Genomic position.}
#' \item{Ref}{Nucleotide(s).}
#' \item{Alt}{Nucleotide(s).}
#' \item{Gene.refGene}{Closest Gene.refGene.}
#' \item{Func.refGene}{Functional consequence of variant.}
#'
#' @export
#' @import dplyr
#' @import stringr
PrepareAnnovarFile <- function(annovar.file) {
  df.annovar.temp <- read.csv(annovar.file, sep = "\t")
  df.annovar.temp$Pos <- df.annovar.temp$Start
  GetClosestGene <- function(x) { # get the closest refGene
    curr.chr <- as.character(x[['Chr']])
    curr.pos <- as.integer(x[['Pos']])
    curr.ref <- toupper(as.character(x[['Ref']]))
    curr.alt <- toupper(as.character(x[['Alt']]))
    curr.gene.refgene <- as.character(x[['Gene.refGene']])
    curr.genedetail.refgene <- as.character(x[['GeneDetail.refGene']])
    curr.func.refgene <- as.character(x[['Func.refGene']])

    # Default data.frame
    df.temp <- data.frame(Chr = curr.chr,
                          Pos = curr.pos,
                          Ref = curr.ref,
                          Alt = curr.alt,
                          Gene.refGene = curr.gene.refgene,
                          Func.refGene = curr.func.refgene,
                          stringsAsFactors = FALSE)

    # Only one Gene.refGene
    curr.gene.refgene <- strsplit(curr.gene.refgene, "\\;")[[1]]
    if (length(curr.gene.refgene) == 1) {
      # Return the original value of Gene.refGene
      return(df.temp)
    }

    # Parse GeneDetail.refGene to return the closest gene
    if (curr.genedetail.refgene == "") { # no "GeneDetail.refGene" value
      # Return the original value of Gene.refGene
      return(df.temp)
    }

    # Return the closest Gene.refGene
    if (str_detect(curr.genedetail.refgene, "dist=") == FALSE) { # no "dist" value
      # Return the original value of Gene.refGene
      return(df.temp)
    }
    curr.genedetail.refgene <- strsplit(curr.genedetail.refgene, "\\;")[[1]]
    curr.genedetail.refgene <- str_remove(curr.genedetail.refgene, "dist=")
    curr.genedetail.refgene[curr.genedetail.refgene == "NONE"] <- "-1"
    refgene.distances <- as.integer(curr.genedetail.refgene)
    closest.refgene.index <- which(refgene.distances == max(refgene.distances))
    if (length(closest.refgene.index) > 1) { # multiple equidistant genes
      # Return the original value of Gene.refGene
      return(df.temp)
    } else {
      df.temp$Gene.refGene <- curr.gene.refgene[closest.refgene.index]
      return(df.temp)
    }
  }

  df.annovar <- apply(df.annovar.temp, MARGIN = 1, FUN = GetClosestGene)
  df.annovar <- bind_rows(df.annovar)
  return(df.annovar)
}
