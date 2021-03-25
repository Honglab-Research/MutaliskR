# Signature identification functions. Adapted from Mutalisk (PMID 29790943).
#
# Author: Andy Jinseok Lee


#' @title Identify contribution weights of underlying mutational signatures
#'
#' @description This function identifies which mutational signatures are underlying a
#' given data.frame and contribution weights of each mutational signature.
#'
#' @param data A data.frame prepared by either PrepareSbsVcfFile, PrepareDbsVcfFile, PrepareIdVcfFile, or PrepareSvVcfFile depending on signature.type.
#' @param reference A data.frame where the first column is the list of mutation types and where 2nd to k-th column headers are the names of the signatures.
#' @param target.signatures A character vector that lists all desired signatures to consider for identification.
#' @param n.cores Number of cores to use (default: 2).
#' @param combn.m An integer value that specifies the number of signatures to consider in each step. 'm' parameter in combn function.
#' @param n.max.signatures An integer value that specifies the maximum number of signatures to consider. It is recommended that n.max.signatures >= initial.exploration.combn.m.
#' @param min.probability A float value that specifies the minimum probability to attribute to a signature.
#' @param zeta.value A float value that is added to the data frequency.
#'
#' @return A data.frame with the following columns:
#' \item{Mutations_Count}{Number of mutations.}
#' \item{Signatures}{Identified mutational signatures separated by comma.}
#' \item{Signatures_Count}{Number of identified mutational signatures.}
#' \item{Signatures_Weights}{Normalized (0 to 1) weights of identified mutational signatures separated by comma.}
#' \item{Mutation_Types}{Mutation types separated by comma.}
#' \item{Mutation_Types_Groups}{Mutation type groups separated by comma.}
#' \item{Observed_Spectrum}{Normalized spectrum (frequency) of observed mutations separated by comma.}
#' \item{Reconstructed_Spectrum}{Normalized spectrum (frequency) of MLE reconstructed mutations separated by comma.}
#' \item{Residual_Spectrum}{Normalized spectrum (frequency) of residual mutations separated by comma.}
#' \item{Cosine_Similarity}{Cosine similarity score between Observed_Spectrum and Reconstructed_Spectrum.}
#' \item{RSS}{Raw residual sum of squares (derived from Residual_Spectrum).}
#' \item{RSS_Normalized}{Normalized residual sum of squares (derived from Residual_Spectrum).}
#' \item{BIC}{Bayesian information criterion of the identified model.}
#'
#' @export
#' @import lsa
#' @import doParallel
#' @import foreach
#' @import parallel
IdentifySignatures <- function(data,
                               reference,
                               target.signatures,
                               n.cores = 2,
                               combn.m = 3,
                               n.max.signatures = 7,
                               min.probability = 0.01,
                               zeta.value = 1e-10) {
  mutation.frequencies <- as.numeric(data[,'Frequency'])
  num.mutations <- sum(as.numeric(mutation.frequencies))
  most.likely.signatures <- c()

  for (i in 1:n.max.signatures) {
    # Combinations of signatures to consider
    signature.combinations <- c()

    if (combn.m > length(target.signatures)) {
      combn.m <- length(target.signatures)
    }

    if (i == 1) {
      # Generate 1, 2, ..., 'combn.m' combinations of signatures
      for (j in 1:combn.m) {
        signature.combinations <- c(signature.combinations, combn(x = target.signatures, m = j, simplify = FALSE))
      }
    } else {
      # Generate combn.m combinations of signatures
      signature.combinations <- combn(x = target.signatures, m = combn.m, simplify = FALSE)
    }

    # Iterate through each combination of signatures by multiprocessing (if available)
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)

    signature.combinations.len <- length(signature.combinations)
    df.mle.results <- foreach (j = 1:signature.combinations.len) %dopar% {
      invisible(most.likely.signatures)
      invisible(signature.combinations)
      invisible(mutation.frequencies)
      invisible(zeta.value)
      invisible(num.mutations)
      invisible(min.probability)

      # Candiate signatures are the most.likely.signatures and the current signature combinations
      candidate.signatures <- c(most.likely.signatures, signature.combinations[j][[1]])
      signature.weights = as.matrix(reference[,candidate.signatures])

      # Perform linear regression by MLE
      mle.results <- EstimateMaximumLikelihood(mutation.frequencies = mutation.frequencies + zeta.value,
                                               signature.weights = signature.weights,
                                               num.mutations = num.mutations,
                                               min.probability = min.probability)

      df.temp <- data.frame(Signatures = paste0(candidate.signatures, collapse = ","),
                            Signatures_Weights = paste0(mle.results$Par, collapse = ","),
                            Cosine_Similarity = mle.results$Cosine_Similarity,
                            RSS = mle.results$RSS,
                            stringsAsFactors = FALSE)
      return(df.temp)
    }

    df.mle.results <- do.call(rbind.data.frame, df.mle.results)
    parallel::stopCluster(cl)

    # Fetch the model with the highest cosine similarity (i.e. lowest rss)
    df.mle.results <- df.mle.results[order(-df.mle.results$Cosine_Similarity),]
    df.mle.results <- df.mle.results[1,]

    # Choose the single most likely signature
    curr.best.signatures <- strsplit(x = df.mle.results$Signatures, split = "\\,")[[1]]
    curr.best.signatures.idx <- which(!(curr.best.signatures %in% most.likely.signatures)) # get indices of signatures not in most.likely.signatures
    curr.best.signatures <- curr.best.signatures[curr.best.signatures.idx] # get signatures not in most.likely.signatres
    curr.best.signatures.weights <- as.numeric(strsplit(x = df.mle.results$Signatures_Weights, split = "\\,")[[1]])
    curr.best.signatures.weights <- curr.best.signatures.weights[curr.best.signatures.idx] # get weights of signatures not in most.likely.signatures
    most.likely.signature <- curr.best.signatures[which.max(curr.best.signatures.weights)]
    most.likely.signatures <- c(most.likely.signatures, most.likely.signature)

    # Log
    if (i == 1) {
      PrintLog(paste0("Selected 1st signature: ", most.likely.signature, "."))
    } else if (i == 2) {
      PrintLog(paste0("Selected 2nd signature: ", most.likely.signature, "."))
    } else if (i == 3) {
      PrintLog(paste0("Selected 3rd signature: ", most.likely.signature, "."))
    } else {
      PrintLog(paste0("Selected ", i, "th signature: ", most.likely.signature, "."))
    }

    # Remove the single most likely signature from target.signatures
    target.signatures <- target.signatures[target.signatures != most.likely.signature]
  }

  # Based on most likely signatures, compute each model
  df.models <- data.frame()
  for (i in 1:length(most.likely.signatures)) {
    curr.model.signatures <- most.likely.signatures[1:i]
    signature.weights = as.matrix(reference[,curr.model.signatures])
    mle.results <- EstimateMaximumLikelihood(mutation.frequencies = mutation.frequencies,
                                             signature.weights = signature.weights,
                                             num.mutations = num.mutations,
                                             min.probability = min.probability)

    # Calculate Bayesian Information Criterion of each model
    # https://en.wikipedia.org/wiki/Bayesian_information_criterion#Gaussian_special_case
    curr.bic.value <- num.mutations * log(mle.results$RSS / num.mutations) + length(curr.model.signatures) * log(num.mutations)

    # Wrap normalized results
    df.temp <- WrapIdentifiedModel(mutations.count = num.mutations,
                                   signatures = curr.model.signatures,
                                   signatures.weights = mle.results$Par_Normalized,
                                   mutation.types = as.character(data[,"Mutation_Type"]),
                                   mutation.types.groups = as.character(data[,"Mutation_Type_Group"]),
                                   observed.spectrum = mutation.frequencies / sum(mutation.frequencies),
                                   reconstructed.spectrum = mle.results$Reconstructed_Spectrum_Normalized,
                                   residual.spectrum = mle.results$Residuals_Normalized,
                                   cosine.similarity = mle.results$Cosine_Similarity,
                                   rss = mle.results$RSS,
                                   rss.normalized = mle.results$RSS_Normalized,
                                   bic = curr.bic.value)

    # Store normalized results
    df.models <- rbind(df.models, df.temp)
  }

  # Sort so that the most parsimonious model is the first row in the data.frame
  df.models <- df.models[order(+df.models$BIC),]

  return(df.models)
}



