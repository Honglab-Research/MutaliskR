# Functions to prepare DBS data for signature identification.
#
# Author: Andy Jinseok Lee


#' @title Preprocess doublet base substitutions
#'
#' @description This function preprocesses a data.frame of DBS.
#'
#' @param df A data.frame with the following columns: Chr, Pos, Ref, Alt.
#' @param bsg BSgenome object.
#' @param analyze.variants A boolean value that indicates whether variant-level signature analysis should be performed.
#'
#' @return A data.frame with the following columns:
#' \item{Chr}{Chromosome name.}
#' \item{Pos}{Genomic position.}
#' \item{Ref}{Reference nucleotide sequence.}
#' \item{Alt}{Alternate nucleotide sequence.}
#' \item{Mutation_Type}{One of the 96 mutation sub types (e.g. AC>CA).}
#'
#' @export
#' @import dplyr
PreprocessDbs <- function(df, bsg, analyze.variants) {
  # Retain doublet base substitutions
  conditions <- ((nchar(df$Ref) == 2) & (nchar(df$Alt) == 2))
  df <- df[conditions,]

  if (nrow(df) == 0) {
    df.temp <- as.data.frame(matrix(, ncol = 6, nrow = 0))
    colnames(df.temp) <- c("Chr", "Pos", "Ref", "Alt",
                           "Mutation_Type", "Mutation_Type_Group")
    return(df.temp)
  }

  df.dbs.signatures <- data.frame(Mutation_Type = DBS.MUTATION.TYPES,
                                  Mutation_Type_Reverse_Complemented = DBS.MUTATION.TYPES.REVERSE.COMPLEMENTED,
                                  Mutation_Type_Group = DBS.MUTATION.TYPES.GROUPS,
                                  stringsAsFactors = FALSE)

  PreprocessDbsHelper <- function(x) {
    curr.chr <- as.character(x[['Chr']])
    curr.pos <- as.integer(x[['Pos']])
    curr.ref <- toupper(as.character(x[['Ref']]))
    curr.alt <- toupper(as.character(x[['Alt']]))
    curr.mutation <- paste0(curr.ref, ">", curr.alt)

    df.reverse.matched <- df.dbs.signatures[df.dbs.signatures$Mutation_Type_Reverse_Complemented == curr.mutation,]
    if (nrow(df.reverse.matched) == 1) { # Take the 'Mutation_Type' value
      curr.mutation <- as.character(df.reverse.matched$Mutation_Type)
    }

    # Fetch the matched row in the reference matrix
    df.matched <- df.dbs.signatures[df.dbs.signatures$Mutation_Type == curr.mutation,]
    if (nrow(df.matched) == 1) {
      curr.mutation.type <- as.character(df.matched$Mutation_Type)
      curr.mutation.type.group <- as.character(df.matched$Mutation_Type_Group)
    } else {
      curr.mutation.type <- ""
      curr.mutation.type.group <- ""
      PrintLog(message = paste0(
        "Could not find the matched 'Mutation_Type' for ",
        curr.chr, ":", curr.pos, ":", curr.ref, ">", curr.alt
        ), type = "WARNING"
      )
    }

    # Prepare a data.frame
    df.temp <- data.frame(Chr = curr.chr,
                          Pos = curr.pos,
                          Ref = curr.ref,
                          Alt = curr.alt,
                          Mutation_Type = curr.mutation.type,
                          Mutation_Type_Group = curr.mutation.type.group,
                          stringsAsFactors = FALSE)

    if (analyze.variants == TRUE) {
      df.temp$Variant_Gene <- as.character(x[["Variant_Gene"]])
      df.temp$Variant_Group <- as.character(x[["Variant_Group"]])
    }
    return(df.temp)
  }

  df.dbs <- apply(df, MARGIN = 1, FUN = PreprocessDbsHelper)
  df.dbs <- bind_rows(df.dbs)
  return(df.dbs)
}

#' @title Prepares DBS data.frame
#'
#' @description
#' This function prepares a DBS data.frame and
#' returns relevant data to run IdentifySignatures function.
#'
#' @param df A data.frame with the following columns: Chr, Pos, Ref, Alt.
#' @param bsg BSgenome object.
#' @param analyze.variants A boolean value that indicates whether variant-level signature analysis should be performed.
#'
#' @return A list with the following items:
#' \item{df.dbs}{A data.frame of DBS data.}
#' \item{df.dbs.frequencies}{A data.frame of DBS mutation tyep frequencies.}
#'
#' @export
#' @import dplyr
PrepareDbsDataFrame <- function(df, bsg, analyze.variants) {
  # Step 1. Preprocess data
  df.dbs <- PreprocessDbs(df = df, bsg = bsg, analyze.variants = analyze.variants)

  # Step 2. Enumerate by mutation type
  df.dbs.frequencies <- data.frame()
  df.dbs.signatures <- data.frame(Mutation_Type = DBS.MUTATION.TYPES,
                                  Mutation_Type_Group = DBS.MUTATION.TYPES.GROUPS,
                                  stringsAsFactors = FALSE)
  for (i in 1:nrow(df.dbs.signatures)) {
    curr.mutation.type <-  df.dbs.signatures[i, "Mutation_Type"]
    curr.mutation.type.group <-  df.dbs.signatures[i, "Mutation_Type_Group"]
    df.temp <- data.frame(Mutation_Type_Group = curr.mutation.type.group,
                          Mutation_Type = curr.mutation.type,
                          Frequency = nrow(df.dbs[df.dbs$Mutation_Type == curr.mutation.type,]),
                          stringsAsFactors = FALSE)
    df.dbs.frequencies <- bind_rows(df.dbs.frequencies, df.temp)
  }
  return(list(df.dbs = df.dbs, df.dbs.frequencies = df.dbs.frequencies))
}

#' @title Prepares DBS VCF file
#'
#' @description
#' This function prepares a DBS VCF file and
#' returns relevant data to run IdentifySignatures function.
#'
#' @param vcf.file VCF file including path.
#' @param bsg BSgenome object.
#'
#' @return A list with the following items:
#' \item{df.dbs}{A data.frame of DBS data.}
#' \item{df.dbs.frequencies}{A data.frame of DBS mutation type frequencies.}
#'
#' @export
#' @import bedr
PrepareDbsVcfFile <- function(vcf.file, bsg) {
  df.temp <- PrepareVcfFile(vcf.file = vcf.file)
  prepared.dbs.data <- PrepareDbsDataFrame(df = df.temp,
                                           bsg = bsg,
                                           analyze.variants = FALSE)
  return(prepared.dbs.data)
}
