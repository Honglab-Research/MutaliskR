# Wrapper functions.
#
# Author: Andy Jinseok Lee


#' @title Wraps an identified model
#'
#' @description This function returns a (wrapped) data.frame of an identified model
#'
#' @param mutations.count Mutations count.
#' @param signatures A character vector of signature names.
#' @param signatures.weights A numeric vector of signature weights.
#' @param mutation.types A character vector of mutation types.
#' @param mutation.types.groups A character vector of mutation type groups.
#' @param observed.spectrum A numeric vector of observed spectrum.
#' @param reconstructed.spectrum A numeric vector of reconstructed spectrum.
#' @param residual.spectrum A numeric vector of residual spectrum.
#' @param cosine.similarity Cosine similiarity score.
#' @param rss Residual sum of squares.
#' @param rss.normalized Normalized residual sum of squares.
#' @param bic Bayesian information criterion.
#'
#' @return A data.frame with the following columns:
#'
#' @export
WrapIdentifiedModel <- function(mutations.count,
                                signatures,
                                signatures.weights,
                                mutation.types,
                                mutation.types.groups,
                                observed.spectrum,
                                reconstructed.spectrum,
                                residual.spectrum,
                                cosine.similarity,
                                rss,
                                rss.normalized,
                                bic) {
  df <- data.frame(Mutations_Count = mutations.count,
                   Signatures = paste0(as.character(signatures), collapse = ","),
                   Signatures_Count = length(signatures),
                   Signatures_Weights = paste0(as.numeric(signatures.weights), collapse = ","),
                   Mutation_Types = paste0(as.character(mutation.types), collapse = ","),
                   Mutation_Types_Groups = paste0(as.character(mutation.types.groups), collapse = ","),
                   Observed_Spectrum = paste0(as.numeric(observed.spectrum), collapse = ","),
                   Reconstructed_Spectrum = paste0(as.numeric(reconstructed.spectrum), collapse = ","),
                   Residual_Spectrum = paste0(as.numeric(residual.spectrum), collapse = ","),
                   Cosine_Similarity = cosine.similarity,
                   RSS = rss,
                   RSS_Normalized = rss.normalized,
                   BIC = bic,
                   stringsAsFactors = FALSE)
  return(df)
}

#' @title Unwraps an identified model
#'
#' @description This function returns an unwrapped list of an identified model.
#'
#' @param identified.model A data.frame returned from
#' \code{\link{IdentifySbsSignatures}}, \code{\link{IdentifyDbsSignatures}}, or \code{\link{IdentifyIndelSignatures}}.
#'
#' @return A data.frame.
#'
#' @export
#' @import stringr
UnwrapIdentifiedModel <- function(identified.model) {
  mutations.count <- identified.model$Mutations_Count
  signatures <- str_split(string = identified.model$Signatures, pattern = "\\,")[[1]]
  signatures.weights <- as.numeric(str_split(string = identified.model$Signatures_Weights, pattern = "\\,")[[1]])
  mutation.types <- str_split(string = identified.model$Mutation_Types, pattern = "\\,")[[1]]
  mutation.types.groups <- str_split(string = identified.model$Mutation_Types_Groups, pattern = "\\,")[[1]]
  observed.spectrum <- as.numeric(str_split(string = identified.model$Observed_Spectrum, pattern = "\\,")[[1]])
  reconstructed.spectrum <- as.numeric(str_split(string = identified.model$Reconstructed_Spectrum, pattern = "\\,")[[1]])
  cosine.similarity <- as.numeric(identified.model$Cosine_Similarity)
  rss <- as.numeric(identified.model$RSS)
  rss.normalized <- as.numeric(identified.model$RSS_Normalized)
  bic <- as.numeric(identified.model$BIC)
  return(list(mutations.count = mutations.count,
              signatures = signatures,
              signatures.weights = signatures.weights,
              mutation.types = mutation.types,
              mutation.types.groups = mutation.types.groups,
              observed.spectrum = observed.spectrum,
              reconstructed.spectrum = reconstructed.spectrum,
              cosine.similarity = cosine.similarity,
              rss = rss,
              rss.normalized = rss.normalized,
              bic = bic))
}
