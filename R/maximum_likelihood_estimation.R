# Maximum likelihood estimation functions. Adapted from Mutalisk (Lee et al., Nucleic Acids Research 2018; PMID 29790943).
#
# Author: Andy Jinseok Lee


#' @title Computes and returns the residual sum of squares
#'
#' @description This function estimates the residual sum of squares.
#'
#' @param par A numeric vector of parameters to optimze.
#' @param mutation.frequencies A numeric vector corresponding to the frequency of each mutation type.
#' @param signature.weights A matrix of signature weights.
#'
#' @return Residual sum of squares.
#'
#' @export
LinearRegressionOptimFn <- function(par, mutation.frequencies, signature.weights) {
  # Spectrum given the current par
  spectrum <- rep(0, length(mutation.frequencies))
  for (i in 1:length(par)) {
    spectrum <- spectrum + (signature.weights[,i] * par[i])
  }

  # Calculate residuals
  residuals <- mutation.frequencies - spectrum

  # Calculate residual sum of squares
  residuals.sum.squared <- sum(residuals^2)

  return(residuals.sum.squared)
}


#' @title Performs maximum likelihood estimation
#'
#' @description This function performs maximum likelihood estimation of a given identification setup.
#'
#' @param mutation.frequencies A numeric vector corresponding to the frequency of each mutation type
#' @param signature.weights A matrix of signature weights.
#' @param num.mutations An integer value that specifies the number of mutations
#' @param min.probability A float value that specifies the minimum probability to attribute to a signature.
#'
#' @return A list with the following items:
#' \item{Cosine_Similarity}{Cosine similarity score between mutation.frequncies and Reconstructed_Spectrum}
#' \item{Par}{A numeric vector of contribution weights of identified signatures.}
#' \item{Par_Normalized}{A normalized (0 to 1) numeric vector of contribution weights of identified signatures.}
#' \item{Reconstructed_Spectrum}{A numeric vector of reconstructed spectrum.}
#' \item{Reconstructed_Spectrum_Normalized}{A normalized (0 to 1) numeric vector of reconstructed spectrum.}
#' \item{Residuals}{A numeric vector of residual spectrum.}
#' \item{Residuals_Normalized}{A normalized (0 to 1) numeric vector of residual spectrum.}
#' \item{RSS}{Residual sum of squares (of Residuals).}
#' \item{RSS_Normalized}{Residual sum of squares (of Residuals_Normalized).}
#'
#' @export
#' @import lsa
EstimateMaximumLikelihood <- function(mutation.frequencies,
                                      signature.weights,
                                      num.mutations,
                                      min.probability = 0.01) {
  # Normalized mutation.frequencies
  mutation.frequencies.normalized <- mutation.frequencies / sum(mutation.frequencies)

  # Get number of signatures
  num.signatures <- ncol(signature.weights)

  # Perform linear regression by maximum likelihood estimation
  result <- optim(par = rep(1, num.signatures),
                  fn = LinearRegressionOptimFn,
                  mutation.frequencies = mutation.frequencies,
                  signature.weights = signature.weights,
                  method = "L-BFGS-B",
                  lower = min.probability * num.mutations,
                  upper = num.mutations)

  # Normalize par
  result$par.normalized <- result$par / sum(result$par)

  # Spectrum given the optimized parameters
  spectrum <- rep(0, length(mutation.frequencies))
  spectrum.normalized <- rep(0, length(mutation.frequencies))
  for (i in 1:length(result$par)) {
    spectrum <- spectrum + (signature.weights[,i] * result$par[i])
    spectrum.normalized <- spectrum.normalized + (signature.weights[,i] * result$par.normalized[i])
  }

  # Calculate residuals
  residuals <- mutation.frequencies - spectrum
  residuals.sum.squared <- sum(residuals^2)

  # Calculate residuals normalized
  residuals.normalized <- mutation.frequencies.normalized - spectrum.normalized
  residuals.normalized.sum.squared <- sum(residuals.normalized^2)

  # Calculate cosine similarity
  cosine.similarity <- lsa::cosine(mutation.frequencies, spectrum)

  return(list(Cosine_Similarity = cosine.similarity,
              Par = result$par,
              Par_Normalized = result$par.normalized,
              Reconstructed_Spectrum = spectrum,
              Reconstructed_Spectrum_Normalized = spectrum.normalized,
              Residuals = residuals,
              Residuals_Normalized = residuals.normalized,
              RSS = residuals.sum.squared,
              RSS_Normalized = residuals.normalized.sum.squared))
}
