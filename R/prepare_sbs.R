# Functions to prepare SBS data for signature identification.
#
# Author: Andy Jinseok Lee


#' @title Preprocess single base substitutions
#'
#' @description This function preprocesses a data.frame of SBS.
#'
#' @param df A data.frame with the following columns: Chr, Pos, Ref, Alt.
#' @param bsg BSgenome object.
#' @param context.length Number of context nucleotides.
#' @param analyze.variants A boolean value that indicates whether variant-level signature analysis should be performed.
#'
#' @return A data.frame with the following columns:
#' \item{Chr}{Chromosome name.}
#' \item{Pos}{Genomic position.}
#' \item{Ref}{Reference nucleotide sequence.}
#' \item{Alt}{Alternate nucleotide sequence.}
#' \item{Type}{One of =the six mutation types (e.g. C>A).}
#' \item{Sub_Type}{One of the 32 trinucleotide subtypes (e.g. ACA).}
#' \item{Somatic_Mutation_Type}{One of the 96 mutation sub types (e.g. A[C>A]A).}
#' \item{Upstream_Seq}{Upstream nucleotide sequence of this SBS.}
#' \item{Downstream_Seq}{Downstream nucleotide sequence of this SBS.}
#'
#' @export
#' @import dplyr
#' @import BSgenome
PreprocessSbs <- function(df,
                          bsg,
                          context.length,
                          analyze.variants) {
  # Retain single base substitutions
  conditions <- ((nchar(df$Ref) == 1) & (df$Ref %in% c("A","T","C","G")) &
                 (nchar(df$Alt) == 1) & (df$Alt %in% c("A","T","C","G")))
  df <- df[conditions,]

  if (nrow(df) == 0) {
    df.temp <- as.data.frame(matrix(, ncol = 10, nrow = 0))
    colnames(df.temp) <- c("Chr",
                           "Pos",
                           "Ref",
                           "Alt",
                           "Mutation_Type_Group",
                           "Mutation_Type_Label",
                           "Mutation_Type")
    return(df.temp)
  }

  offset <- (context.length - 1) / 2 # e.g. 1 if tri-nucleotide, 2 if penta-nucleotide

  PreprocessSbsHelper <- function(x) {
    curr.chr <- as.character(x[['Chr']])
    curr.pos <- as.integer(x[['Pos']])
    curr.ref <- toupper(as.character(x[['Ref']]))
    curr.alt <- toupper(as.character(x[['Alt']]))

    # Convert to pyrimidines
    if (curr.ref %in% c("A","G")) {
      # Reference sequence is purine (A, G; take complement)
      curr.ref <- GetComplementBase(nucleotide = curr.ref)
      curr.alt <- GetComplementBase(nucleotide = curr.alt)
      is.ref.purine <- TRUE
      downstream.start <- curr.pos - offset
      downstream.end <- curr.pos - 1
      upstream.start <- curr.pos + 1
      upstream.end <- curr.pos + offset
    } else {
      # Reference sequence is pyrimidine (T, C)
      is.ref.purine <- FALSE
      downstream.start <- curr.pos + 1
      downstream.end <- curr.pos + offset
      upstream.start <- curr.pos - offset
      upstream.end <- curr.pos - 1
    }

    # Prepare a data.frame
    df.temp <- data.frame(Chr = curr.chr,
                          Pos = curr.pos,
                          Ref = curr.ref,
                          Alt = curr.alt,
                          Is_Ref_Purine = is.ref.purine,
                          Upstream_Start = upstream.start,
                          Upstream_End = upstream.end,
                          Downstream_Start = downstream.start,
                          Downstream_End = downstream.end,
                          stringsAsFactors = FALSE)
    if (analyze.variants == TRUE) {
      df.temp$Variant_Gene <- as.character(x[["Variant_Gene"]])
      df.temp$Variant_Group <- as.character(x[["Variant_Group"]])
    }
    return(df.temp)
  }

  df.sbs <- apply(df, MARGIN = 1, FUN = PreprocessSbsHelper)
  df.sbs <- bind_rows(df.sbs)

  # Get upstream and downstream sequences
  df.sbs$Upstream_Seq <- as.character(
    BSgenome::getSeq(bsg,
                     names = df.sbs$Chr,
                     start = df.sbs$Upstream_Start,
                     end = df.sbs$Upstream_End)
  )
  df.sbs$Downstream_Seq <- as.character(
    BSgenome::getSeq(bsg,
                     names = df.sbs$Chr,
                     start = df.sbs$Downstream_Start,
                     end = df.sbs$Downstream_End)
  )

  # Convert purines to pyrimidines
  if (context.length == 3) {
    df.sbs[df.sbs$Is_Ref_Purine, "Upstream_Seq"] <- sapply(df.sbs[df.sbs$Is_Ref_Purine, "Upstream_Seq"], FUN = GetComplementBase)
    df.sbs[df.sbs$Is_Ref_Purine, "Downstream_Seq"] <- sapply(df.sbs[df.sbs$Is_Ref_Purine, "Downstream_Seq"], FUN = GetComplementBase)
  } else {
    df.sbs[df.sbs$Is_Ref_Purine, "Upstream_Seq"] <- sapply(df.sbs[df.sbs$Is_Ref_Purine, "Upstream_Seq"], FUN = GetComplementBases)
    df.sbs[df.sbs$Is_Ref_Purine, "Downstream_Seq"] <- sapply(df.sbs[df.sbs$Is_Ref_Purine, "Downstream_Seq"], FUN = GetComplementBases)
  }

  # Create labels
  df.sbs$Mutation_Type_Group <- paste0(df.sbs$Ref, ">", df.sbs$Alt) # e.g. C>A
  df.sbs$Mutation_Type_Label <- paste0(df.sbs$Upstream_Seq, df.sbs$Ref, df.sbs$Downstream_Seq) # e.g. ACA
  df.sbs$Mutation_Type <- paste0(df.sbs$Upstream_Seq,
                                 "[",df.sbs$Ref, ">", df.sbs$Alt, "]",
                                 df.sbs$Downstream_Seq) # e.g. A[C>A]A
  return(df.sbs)
}

#' @title Prepares SBS data.frame
#'
#' @description
#' This function prepares a SBS data.frame and
#' returns relevant data to run IdentifySignatures function.
#'
#' @param df A data.frame with the following columns: Chr, Pos, Ref, Alt.
#' @param bsg BSgenome object.
#' @param reference A data.frame with the following columns: Mutation_Type, and names of SBS signatures (default: a data.frame returned from \code{\link{GetPcawgSbsSignaturesData}}).
#' @param context.length Number of context nucleotides.
#' @param analyze.variants A boolean value that indicates whether variant-level signature analysis should be performed (default: FALSE).
#'
#' @return A list with the following items:
#' \item{df.sbs}{A data.frame of SBS data.}
#' \item{df.sbs.frequencies}{A data.frame of SBS mutation tyep frequencies.}
#'
#' @export
#' @import dplyr
PrepareSbsDataFrame <- function(df,
                                bsg,
                                reference,
                                context.length,
                                analyze.variants) {
  # Step 1. Preprocess data
  df.sbs <- PreprocessSbs(df = df,
                          bsg = bsg,
                          context.length = context.length,
                          analyze.variants = analyze.variants)

  # Step 2. Enumerate count by mutation type
  if (context.length == 3) {
    df.sbs.signatures <- data.frame(Mutation_Type = SBS.MUTATION.TYPES,
                                    Mutation_Type_Group = SBS.MUTATION.TYPES.GROUPS,
                                    stringsAsFactors = FALSE)
  } else {
    df.sbs.signatures <- data.frame(Mutation_Type = reference$Mutation_Type,
                                    Mutation_Type_Group = "",
                                    stringsAsFactors = FALSE)
  }
  df.sbs.frequencies <- data.frame()
  for (i in 1:nrow(df.sbs.signatures)) {
    curr.mutation.type <- as.character(df.sbs.signatures[i, "Mutation_Type"])
    curr.mutation.type.group <- as.character(df.sbs.signatures[i, "Mutation_Type_Group"])
    df.temp <- data.frame(Mutation_Type = curr.mutation.type,
                          Mutation_Type_Group = curr.mutation.type.group,
                          Frequency = nrow(df.sbs[df.sbs$Mutation_Type == curr.mutation.type,]),
                          stringsAsFactors = FALSE)
    df.sbs.frequencies <- bind_rows(df.sbs.frequencies, df.temp)
  }
  return(list(df.sbs = df.sbs, df.sbs.frequencies = df.sbs.frequencies))
}

#' @title Prepares SBS VCF file
#'
#' @description
#' This function prepares a SBS VCF file and
#' returns relevant data to run IdentifySignatures function.
#'
#' @param vcf.file VCF file including path.
#' @param bsg BSgenome object.
#' @param reference A data.frame with the following columns: Mutation_Type, and names of SBS signatures (default: a data.frame returned from \code{\link{GetPcawgSbsSignaturesData}}).
#' @param context.length Number of context nucleotides.
#'
#' @return A list with the following items:
#' \item{df.sbs}{A data.frame of SBS data.}
#' \item{df.sbs.frequencies}{A data.frame of SBS mutation type frequencies.}
#'
#' @export
#' @import bedr
PrepareSbsVcfFile <- function(vcf.file, bsg, reference, context.length) {
  df.temp <- PrepareVcfFile(vcf.file = vcf.file)
  prepared.sbs.data <- PrepareSbsDataFrame(df = df.temp,
                                           bsg = bsg,
                                           reference = reference,
                                           context.length = context.length,
                                           analyze.variants = FALSE)
  return(prepared.sbs.data)
}
