# Variant signature assignment functions. Adapted from Letouze et al. (PMID 29101368).
#
# Author: Andy Jinseok Lee


#' @title Identifies mutational signatures associated with each variant
#'
#' @description This function identifies mutational signatures associated with each variant.
#'
#' @param prepared.data A data.frame returned from \code{\link{PrepareSbsDataFrame}}.
#' @param identified.model A data.frame returned from \code{\link{IdentifySignatures}}.
#' @param reference A data.frame with the following columns: Mutation_Type, and names of SBS signatures.
#'
#' @return A data.frame with the following columns:
#' \item{Chr}{Chromosome.}
#' \item{Pos}{Genomic position.}
#' \item{Ref}{Reference allele.}
#' \item{Alt}{Alternate allele.}
#' \item{Mutation_Type}{Mutation type.}
#' \item{Variant_Gene}{Variant gene.}
#' \item{Variant_Group}{Variant group.}
#' \item{Signature}{Signature.}
#' \item{Probability}{Probability the mutation is the result of the signature.}
#'
#' @export
#' @import dplyr
IdentifyVariantSignatures <- function(prepared.data, identified.model, reference) {
  df <- prepared.data[[1]]
  df.frequencies <- prepared.data[[2]]
  unwrapped.identified.model <- UnwrapIdentifiedModel(identified.model = identified.model)

  ComputeProbabilityByMutationType <- function(x) {
    # x =
    curr.mutation.type <- as.character(x[['Mutation_Type']])
    n.c <- nrow(df.frequencies[df.frequencies$Mutation_Type == curr.mutation.type,])

    # Calculate probability of current mutation due to each signature identified
    df.curr.mutation.type.sigs.probs <- data.frame()
    sigma.value <- 0
    for (i in 1:length(unwrapped.identified.model$signatures)) {
      curr.signature <- unwrapped.identified.model$signatures[i]
      w.j <- unwrapped.identified.model$signatures.weights[i]
      r.c.j <- as.numeric(reference[reference$Mutation_Type == curr.mutation.type, curr.signature])
      curr.signature.prob <- n.c * w.j * r.c.j
      df.temp <- data.frame(Mutation_Type = curr.mutation.type,
                            Signature = curr.signature,
                            Probability = curr.signature.prob,
                            stringsAsFactors = FALSE)
      df.curr.mutation.type.sigs.probs <- bind_rows(df.curr.mutation.type.sigs.probs, df.temp)
      sigma.value <- sigma.value + (n.c * w.j * r.c.j)
    }

    df.curr.mutation.type.sigs.probs$Probability <- df.curr.mutation.type.sigs.probs$Probability / sigma.value
    return(df.curr.mutation.type.sigs.probs)
  }

  df.mut.types.sigs.probs <- apply(df.frequencies, MARGIN = 1, FUN = ComputeProbabilityByMutationType)
  df.mut.types.sigs.probs <- bind_rows(df.mut.types.sigs.probs)

  IdentifyVariantSignaturesHelper <- function(x) {
    curr.chr <- as.character(x[['Chr']])
    curr.pos <- as.character(x[['Pos']])
    curr.ref <- as.character(x[['Ref']])
    curr.alt <- as.character(x[['Alt']])
    curr.mutation.type <- as.character(x[['Mutation_Type']])
    curr.gene <- as.character(x[["Variant_Gene"]])
    curr.group <- as.character(x[["Variant_Group"]])

    df.temp <- df.mut.types.sigs.probs[df.mut.types.sigs.probs$Mutation_Type == curr.mutation.type,]
    df.temp$Chr <- curr.chr
    df.temp$Pos <- curr.pos
    df.temp$Ref <- curr.ref
    df.temp$Alt <- curr.alt
    df.temp$Variant_Gene <- curr.gene
    df.temp$Variant_Group <- curr.group
    return(df.temp)
  }

  df.sigs.probs <- apply(df, MARGIN = 1, FUN = IdentifyVariantSignaturesHelper)
  df.sigs.probs <- bind_rows(df.sigs.probs)
  df.sigs.probs <- df.sigs.probs[,c("Chr",
                                    "Pos",
                                    "Ref",
                                    "Alt",
                                    "Mutation_Type",
                                    "Variant_Gene",
                                    "Variant_Group",
                                    "Signature",
                                    "Probability")]
  return(df.sigs.probs)
}
