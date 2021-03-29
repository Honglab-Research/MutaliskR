# Main functions.
#
# Author: Andy Jinseok Lee


#' @title Identify single-base substitution mutational signatures
#'
#' @description This function identifies SBS mutational signatures.
#'
#' @param input Either VCF file path or data.frame. If the input is a data.frame, it must include the following columns: Chr, Pos, Ref, Alt.
#' @param bsg BSgenome object.
#' @param sample.id Sample ID that will be used to name output files (default: 'Sample').
#' @param reference A data.frame with the following columns: Mutation_Type, and names of SBS signatures (default: a data.frame returned from \code{\link{GetPcawgSbsSignaturesData(sequencing.type = "WGS")}}).
#' @param target.signatures Signatures to be considered for identification (default: an array returned from \code.{\link{GetPcawgSbsSignaturesNames(sequencing.type = "WGS")}}).
#' @param plot.theme A data.frame returned from \code{\link{GetSbsSignaturesPlotTheme}}.
#' @param analyze.variants.column.gene Name of column in the data.frame corresponding to the gene name or ID (e.g. "Gene.refGene" if using ANNOVAR data).
#' @param analyze.variants.column.group Name of column in the data.frame corresponding to the variant group for plotting purposes (e.g. "Func.refGene" if using ANNOVAR data).
#' @param analyze.variants A boolean value that indicates whether variant-level signature analysis should be performed (default: FALSE).
#' @param context.length Number of context nucleotides (default: 3; tri-nucleotide context).
#' @param n.cores Number of cores to use.
#' @param combn.m Number of signatures to consider in each step. 'm' parameter in combn function (default: 3).
#' @param n.max.signatures Maximum number of signatures to identify. Recommended: n.max.signatures >= initial.exploration.combn.m (default: 7).
#' @param min.probability Minimum probability to attribute to a signature (default: 0.01).
#' @param zeta.value A float value that is added to the data frequency (default: 1e-10).
#' @param save Save resulting files if TRUE, otherwise do not save (default: TRUE).
#' @param save.dir Save directory path (default: NULL).
#'
#' @return A list with the following elements:
#' @return results: a data.frame with the following columns:
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
#' @return results.variants: a data.frame with the following columns:
#'
#' @export
#' @import dplyr
IdentifySbsSignatures <- function(input,
                                  bsg,
                                  sample.id = "Sample",
                                  reference = GetPcawgSbsSignaturesData(sequencing.type = "WGS"),
                                  target.signatures = GetPcawgSbsSignaturesNames(sequencing.type = "WGS"),
                                  plot.theme = GetSbsSignaturesPlotTheme(),
                                  analyze.variants.column.gene,
                                  analyze.variants.column.group,
                                  analyze.variants = FALSE,
                                  context.length = 3,
                                  n.cores = 2,
                                  combn.m = 3,
                                  n.max.signatures = 7,
                                  min.probability = 0.01,
                                  zeta.value = 1e-10,
                                  save = TRUE,
                                  save.dir = NULL) {
  # Step 1. Check integrity of inputs
  PrintLog("Step 1/6 - Checking integrity of input parameters.")
  if (is.data.frame(input) == TRUE) { # input is data.frame
    input.type <- "data.frame"
  } else {                            # input is VCF file
    if (file.exists(input) == TRUE) {
      input.type <- "file"
    } else {
      PrintLog("The parameter 'input' appears to be a file but does not exist.", type = "ERROR")
      return()
    }
    if (analyze.variants == TRUE) {
      PrintLog("The parameter 'analyze.variants' cannot be TRUE if the input is a VCF file", type = "ERROR")
      return()
    }
  }
  if (save == TRUE) {
    if (is.null(save.dir) == TRUE) {
      PrintLog("The parameter 'save.dir' must not be null if the parameter 'save' is TRUE.", type = "ERROR")
      return()
    } else {
      if (identical(strsplit(x = save.dir, split = "")[[1]], "/") == FALSE) {
        save.dir <- paste0(save.dir, "/")
      }
      if (dir.exists(save.dir) == FALSE) {
        dir.create(path = save.dir, recursive = TRUE, showWarnings = FALSE)
      }
    }
  }
  if (analyze.variants == TRUE) {
    if (is.null(analyze.variants.column.gene)) {
      PrintLog("The parameter 'analyze.variants.column.gene' must not be null if the parameter 'analyze.variants' is TRUE.", type = "ERROR")
      return()
    }
    if (is.null(analyze.variants.column.group)) {
      PrintLog("The parameter 'analyze.variants.column.group' must not be null if the parameter 'analyze.variants' is TRUE.", type = "ERROR")
      return()
    }
  }
  if (is.numeric(n.cores) == FALSE) {
    PrintLog("The parameter 'n.cores' must be a number.", type = "ERROR")
    return()
  } else {
    if (n.cores <= 0) {
      PrintLog("The parameter 'n.cores' must be a positive number", type = "ERROR")
      return()
    }
  }
  if (is.numeric(context.length) == FALSE) {
    PrintLog("The parameter 'context.length' must be a number", type = "ERROR")
    return()
  } else {
    if (context.length < 3 || (context.length %% 2) == 0) {
      PrintLog("The parameter 'context.length' must be an odd number equal to or greater than 3.", type = "ERROR")
      return()
    }
  }
  if (is.numeric(combn.m) == FALSE) {
    PrintLog("The parameter 'combn.m' must be a number.", type = "ERROR")
    return()
  }
  if (is.numeric(n.max.signatures) == FALSE) {
    PrintLog("The parameter 'n.max.signatures' must be a number.", type = "ERROR")
    return()
  }
  if (is.numeric(min.probability) == FALSE) {
    PrintLog("The parameter 'min.probability' must be a number.", type = "ERROR")
    return()
  }
  if (is.numeric(zeta.value) == FALSE) {
    PrintLog("The parameter 'zeta.value' must be a number.", type = "ERROR")
    return()
  }

  # Step 2. Preprocess input variants
  PrintLog("Step 2/6 - Preprocessing input variants.")
  if (analyze.variants) {
    input$Variant_Gene <- input[[`analyze.variants.column.gene`]]
    input$Variant_Group <- input[[`analyze.variants.column.group`]]
  }
  if (input.type == "data.frame") {
    prepared.sbs.data <- PrepareSbsDataFrame(
      df = input,
      bsg = bsg,
      reference = reference,
      context.length = context.length,
      analyze.variants = analyze.variants
    )
  }
  if (input.type == "file") {
    prepared.sbs.data <- PrepareSbsVcfFile(
      vcf.file = input,
      bsg = bsg,
      reference = reference,
      context.length = context.length
    )
  }

  # Step 3. Write prepared SBS data to files
  if (save == TRUE) {
    PrintLog("Step 3/6 - Writing preprocessed SBS data to files.")
    write.table(x = prepared.sbs.data$df.sbs,
                file = paste0(save.dir, sample.id, "_SBS_Data.tsv"),
                sep = "\t", row.names = FALSE)
    write.table(x = prepared.sbs.data$df.sbs.frequencies,
                file = paste0(save.dir, sample.id, "_SBS_Frequency.tsv"),
                sep = "\t", row.names = FALSE)

    # Create sub-directory for the i-th signature model
    for (i in 1:n.max.signatures) {
      dir.create(paste0(save.dir, sample.id, "_", i, "_Signature_Model/"))
    }
  } else {
    PrintLog("Step 3/6 - SKIP writing preprocessed SBS data to files.")
  }

  # Step 4. Identify signatures
  PrintLog("Step 4/6 - Identifying SBS mutational signatures.")
  results <- IdentifySignatures(data = prepared.sbs.data$df.sbs.frequencies,
                                reference = reference,
                                target.signatures = target.signatures,
                                n.cores = n.cores,
                                combn.m = combn.m,
                                n.max.signatures = n.max.signatures,
                                min.probability = min.probability,
                                zeta.value = zeta.value)
  results$Sample_ID <- sample.id

  # Step 5. Identify signatures per variant
  results.variants <- data.frame()
  if (analyze.variants == TRUE) {
    PrintLog("Step 5/6 - Identifying SBS mutational signatures for each variant.")
    for (i in 1:nrow(results)) {
      signatures.count <- as.integer(results[i, 'Signatures_Count'])
      df.curr.model.variant.signatures <- IdentifyVariantSignatures(prepared.data = prepared.sbs.data,
                                                                    identified.model = results[i,],
                                                                    reference = reference)
      df.curr.model.variant.signatures$Signatures_Count <- signatures.count
      results.variants <- bind_rows(results.variants, df.curr.model.variant.signatures)
    }
  } else {
    PrintLog("Step 5/6 - SKIP identification of SBS mutational signatures for each variant.")
  }

  # Step 6. Save data
  if (save == TRUE) {
    PrintLog("Step 6/6 - Saving data.")
    write.table(x = results,
                file = paste0(save.dir, sample.id, "_All_Models.tsv"),
                sep = "\t", row.names = FALSE)

    if (analyze.variants == TRUE) {
      write.table(x = results.variants,
                  file = paste0(save.dir, sample.id, "_All_Models_Variant_Signatures.tsv"),
                  sep = "\t", row.names = FALSE)
    }

    # All models
    for (i in 1:nrow(results)) {
      # Signatures
      if (context.length == 3) {
        signatures.count <- as.integer(results[i, 'Signatures_Count'])
        curr.model.plots <- PlotIdentifiedModel(identified.model = results[i,],
                                                reference = reference,
                                                plot.theme = plot.theme)
        ggsave(plot = curr.model.plots$plot.merged,
               filename = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                                 sample.id, "_", signatures.count, "_Signature_Model.svg"),
               width = 16, height = 10, dpi = 600)

        curr.model.weight.plot <- PlotIdentifiedModelSignatureWeights(identified.model = results[i,],
                                                                      df.signatures.colors = GetSbsSignaturesColors(),
                                                                      plot.theme = plot.theme)
        ggsave(plot = curr.model.weight.plot,
               filename = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                                 sample.id, "_", signatures.count, "_Signature_Model_Weights.svg"),
               width = 16, height = 10, dpi = 600)

        # Signatures by variant
        if (analyze.variants == TRUE) {
          df.curr.variant.signatures.model <- results.variants[results.variants$Signatures_Count == signatures.count,]
          write.table(x = df.curr.variant.signatures.model,
                      file = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                                    sample.id, "_", signatures.count, "_Signature_Model_Variant_Signatures.tsv"),
                      sep = "\t", row.names = FALSE)
          plots.variants <- PlotIdentifiedModelVariantSignatures(df.sigs.probs = df.curr.variant.signatures.model,
                                                                 df.signatures.colors = GetSbsSignaturesColors())
          # Signature by variant plots
          if (length(plots.variants$signature.group.plots) > 0) { # signature group
            for(j in 1:length(plots.variants$signature.group.plots)) {
              ggsave(plot = plots.variants$signature.group.plots[[j]],
                     filename = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                                       sample.id, "_", signatures.count, "_Signature_Model_Variant_Signature_Groups_",
                                       plots.variants$signature.group.plots.names[j], ".svg"),
                     width = 16, height = 9, dpi = 600)
            }
          }
          if (length(plots.variants$signature.plots) > 0) { # signature
            for(j in 1:length(plots.variants$signature.plots)) {
              ggsave(plot = plots.variants$signature.plots[[j]],
                     filename = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                                       sample.id, "_", signatures.count, "_Signature_Model_Variant_Signatures_",
                                       plots.variants$signature.plots.names[j], ".svg"),
                     width = 16, height = 9, dpi = 600)
            }
          }
        }
      } else {
        PrintLog(message = "SBS plot can only be created for tri-nucleotide context", type = "WARNING")
      }
    }

    # Most parsimonious model
    # Signatures
    if (context.length == 3) {
      plots.models <- PlotIdentifiedModel(identified.model = results[1,],
                                          reference = reference,
                                          plot.theme = plot.theme)
      ggsave(plot = plots.models$plot.merged,
             filename = paste0(save.dir, sample.id, "_Best_Model.svg"),
             width = 16, height = 10, dpi = 600)
    }

    write.table(x = results[1,],
                file = paste0(save.dir, sample.id, "_Best_Model.tsv"),
                sep = "\t", row.names = FALSE)

    if (context.length == 3) {
      plot.weights <- PlotIdentifiedModelSignatureWeights(identified.model = results[1,],
                                                          df.signatures.colors = GetSbsSignaturesColors(),
                                                          plot.theme = plot.theme)
      ggsave(plot = plot.weights,
             filename = paste0(save.dir, sample.id, "_Best_Model_Weights.svg"),
             width = 16, height = 10, dpi = 600)
    }

    # Signatures by variant
    if (analyze.variants == TRUE) {
      best.model.signatures.count <- as.integer(results[1, 'Signatures_Count'])
      best.model.variant.signatures <- results.variants[results.variants$Signatures_Count == best.model.signatures.count,]

      if (context.length == 3) {
        plots.variants <- PlotIdentifiedModelVariantSignatures(df.sigs.probs = best.model.variant.signatures,
                                                               df.signatures.colors = GetSbsSignaturesColors())
        if (length(plots.variants$signature.group.plots) > 0) { # signature group
          for(j in 1:length(plots.variants$signature.group.plots)) {
            ggsave(plot = plots.variants$signature.group.plots[[j]],
                   filename = paste0(save.dir, sample.id, "_Best_Model_Variant_Signature_Groups_",
                                     plots.variants$signature.group.plots.names[j], ".svg"),
                   width = 16, height = 9, dpi = 600)
          }
        }
        if (length(plots.variants$signature.plots) > 0) { # signature
          for(j in 1:length(plots.variants$signature.plots)) {
            ggsave(plot = plots.variants$signature.plots[[j]],
                   filename = paste0(save.dir, sample.id, "_Best_Model_Variant_Signatures_",
                                     plots.variants$signature.plots.names[j], ".svg"),
                   width = 16, height = 9, dpi = 600)
          }
        }
      }

      write.table(x = best.model.variant.signatures,
                  file = paste0(save.dir, sample.id, "_Best_Model_Variant_Signatures.tsv"),
                  sep = "\t", row.names = FALSE)
    }
  } else {
    PrintLog("Step 6/6 - SKIP saving data.")
  }

  return(list(results = results,
              results.variants = results.variants))
}

#' @title Identify doublet base subsitution mutational signatures
#'
#' @description This function identifies DBS mutational signatures.
#'
#' @param input Either VCF file path or data.frame. If the input is a data.frame, it must include the following columns: Chr, Pos, Ref, Alt.
#' @param bsg BSgenome object.
#' @param sample.id Sample ID that will be used to name output files (default: 'Sample').
#' @param reference A data.frame with the following columns: Mutation_Type, and names of DBS signatures
#' (default: a data.frame returned from \code{\link{GetPcawgDbsSignaturesData}}).
#' @param target.signatures Signatures to be considered for identification (default: an array returned from \code{\link{GetPcawgDbsSignaturesNames}}).
#' @param plot.theme A data.frame returned from \code{\link{GetDbsSignaturesPlotTheme}}.
#' @param analyze.variants.column.gene Name of column in the data.frame corresponding to the gene name or ID (e.g. "Gene.refGene" if using ANNOVAR data).
#' @param analyze.variants.column.group Name of column in the data.frame corresponding to the variant group for plotting purposes (e.g. "Func.refGene" if using ANNOVAR data).
#' @param n.cores Number of cores to use.
#' @param combn.m Number of signatures to consider in each step. 'm' parameter in combn function (default: 3).
#' @param n.max.signatures Maximum number of signatures to identify. Recommended: n.max.signatures >= initial.exploration.combn.m (default: 7).
#' @param min.probability Minimum probability to attribute to a signature (default: 0.01).
#' @param zeta.value A float value that is added to the data frequency (default: 1e-10).
#' @param save Save resulting files if TRUE, otherwise do not save (default: TRUE).
#' @param save.dir Save directory path (default: NULL).
#'
#' @return A list with the following elements:
#' @return results: a data.frame with the following columns:
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
#' @return results.variants: a data.frame with the following columns:
#'
#' @export
#' @import dplyr
IdentifyDbsSignatures <- function(input,
                                  bsg,
                                  sample.id = "Sample",
                                  reference = GetPcawgDbsSignaturesData(),
                                  target.signatures = GetPcawgDbsSignaturesNames(),
                                  plot.theme = GetDbsSignaturesPlotTheme(),
                                  analyze.variants.column.gene,
                                  analyze.variants.column.group,
                                  analyze.variants = FALSE,
                                  n.cores = 2,
                                  combn.m = 3,
                                  n.max.signatures = 7,
                                  min.probability = 0.01,
                                  zeta.value = 1e-10,
                                  save = TRUE,
                                  save.dir = NULL) {
  # Step 1. Check integrity of inputs
  PrintLog("Step 1/6 - Checking integrity of input parameters.")
  if (is.data.frame(input) == TRUE) { # input is data.frame
    input.type <- "data.frame"
  } else {
    if (file.exists(input) == TRUE) { # input is VCF file
      input.type <- "file"
    } else {
      PrintLog("The parameter 'input' appears to be a file but does not exist.", type = "ERROR")
      return()
    }
    if (analyze.variants == TRUE) {
      PrintLog("The parameter 'analyze.variants' cannot be TRUE if the input is a VCF file", type = "ERROR")
      return()
    }
  }
  if (save == TRUE) {
    if (is.null(save.dir) == TRUE) {
      PrintLog("The parameter 'save.dir' must be set if the parameter 'save' is TRUE.", type = "ERROR")
      return()
    } else {
      if (identical(strsplit(x = save.dir, split = "")[[1]], "/") == FALSE) {
        save.dir <- paste0(save.dir, "/")
      }
      if (dir.exists(save.dir) == FALSE) {
        dir.create(path = save.dir, recursive = TRUE, showWarnings = FALSE)
      }
    }
  }
  if (analyze.variants == TRUE) {
    if (is.null(analyze.variants.column.gene)) {
      PrintLog("The parameter 'analyze.variants.column.gene' must not be null if the parameter 'analyze.variants' is TRUE.", type = "ERROR")
      return()
    }
    if (is.null(analyze.variants.column.group)) {
      PrintLog("The parameter 'analyze.variants.column.group' must not be null if the parameter 'analyze.variants' is TRUE.", type = "ERROR")
      return()
    }
  }
  if (is.numeric(n.cores) == FALSE) {
    PrintLog("The parameter 'n.cores' must be a number.", type = "ERROR")
    return()
  } else {
    if (n.cores <= 0) {
      PrintLog("The parameter 'n.cores' must be a positive number", type = "ERROR")
      return()
    }
  }
  if (is.numeric(combn.m) == FALSE) {
    PrintLog("The parameter 'combn.m' must be a number.", type = "ERROR")
    return()
  }
  if (is.numeric(n.max.signatures) == FALSE) {
    PrintLog("The parameter 'n.max.signatures' must be a number.", type = "ERROR")
    return()
  }
  if (is.numeric(min.probability) == FALSE) {
    PrintLog("The parameter 'min.probability' must be a number.", type = "ERROR")
    return()
  }
  if (is.numeric(zeta.value) == FALSE) {
    PrintLog("The parameter 'zeta.value' must be a number.", type = "ERROR")
    return()
  }

  # Step 2. Preprocess input data
  PrintLog("Step 2/6 - Preprocessing input variants.")
  if (analyze.variants) {
    input$Variant_Gene <- input[[`analyze.variants.column.gene`]]
    input$Variant_Group <- input[[`analyze.variants.column.group`]]
  }
  if (input.type == "data.frame") {
    prepared.dbs.data <- PrepareDbsDataFrame(
      df = input,
      bsg = bsg,
      analyze.variants = analyze.variants
    )
  }
  if (input.type == "file") {
    prepared.dbs.data <- PrepareDbsVcfFile(
      vcf.file = input,
      bsg = bsg
    )
  }

  # Step 3. Write prepared DBS data to files
  if (save == TRUE) {
    PrintLog("Step 3/6 - Writing preprocessed DBS data to files.")
    write.table(x = prepared.dbs.data$df.dbs,
                file = paste0(save.dir, sample.id, "_DBS_Data.tsv"),
                sep = "\t", row.names = FALSE)
    write.table(x = prepared.dbs.data$df.dbs.frequencies,
                file = paste0(save.dir, sample.id, "_DBS_Frequency.tsv"),
                sep = "\t", row.names = FALSE)

    # Create sub-directory for the i-th signature model
    for (i in 1:n.max.signatures) {
      dir.create(paste0(save.dir, sample.id, "_", i, "_Signature_Model/"))
    }
  } else {
    PrintLog("Step 3/6 - SKIP writing preprocessed DBS data to files.")
  }

  # Step 4. Identify signatures
  PrintLog("Step 4/6 - Identifying DBS mutational signatures.")
  results <- IdentifySignatures(data = prepared.dbs.data$df.dbs.frequencies,
                                reference = reference,
                                target.signatures = target.signatures,
                                n.cores = n.cores,
                                combn.m = combn.m,
                                n.max.signatures = n.max.signatures,
                                min.probability = min.probability,
                                zeta.value = zeta.value)
  results$Sample_ID <- sample.id

  # Step 5. Identify signatures per variant
  results.variants <- data.frame()
  if (analyze.variants == TRUE) {
    PrintLog("Step 5/6 - Identifying DBS mutational signatures for each variant.")
    for (i in 1:nrow(results)) {
      signatures.count <- as.integer(results[i, 'Signatures_Count'])
      df.curr.model.variant.signatures <- IdentifyVariantSignatures(prepared.data = prepared.dbs.data,
                                                                    identified.model = results[i,],
                                                                    reference = reference)
      df.curr.model.variant.signatures$Signatures_Count <- signatures.count
      results.variants <- bind_rows(results.variants, df.curr.model.variant.signatures)
    }
  } else {
    PrintLog("Step 5/6 - SKIP identification of DBS mutational signatures for each variant.")
  }

  # Step 6. Save data
  if (save == TRUE) {
    PrintLog("Step 6/6 - Saving data.")
    write.table(x = results,
                file = paste0(save.dir, sample.id, "_All_Models.tsv"),
                sep = "\t", row.names = FALSE)

    if (analyze.variants == TRUE) {
      write.table(x = results.variants,
                  file = paste0(save.dir, sample.id, "_All_Models_Variant_Signatures.tsv"),
                  sep = "\t", row.names = FALSE)
    }

    # All models
    for (i in 1:nrow(results)) {
      # Signatures
      signatures.count <- as.integer(results[i, 'Signatures_Count'])
      curr.model.plots <- PlotIdentifiedModel(identified.model = results[i,],
                                              reference = reference,
                                              plot.theme = plot.theme)
      ggsave(plot = curr.model.plots$plot.merged,
             filename = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                               sample.id, "_", signatures.count, "_Signature_Model.svg"),
             width = 16, height = 10, dpi = 600)

      curr.model.weight.plot <- PlotIdentifiedModelSignatureWeights(identified.model = results[i,],
                                                                    df.signatures.colors = GetDbsSignaturesColors(),
                                                                    plot.theme = plot.theme)
      ggsave(plot = curr.model.weight.plot,
             filename = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                               sample.id, "_", signatures.count, "_Signature_Model_Weights.svg"),
             width = 16, height = 10, dpi = 600)

      # Signatures by variant
      if (analyze.variants == TRUE) {
        df.curr.variant.signatures.model <- results.variants[results.variants$Signatures_Count == signatures.count,]
        write.table(x = df.curr.variant.signatures.model,
                    file = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                                  sample.id, "_", signatures.count, "_Signature_Model_Variant_Signatures.tsv"),
                    sep = "\t", row.names = FALSE)
        plots.variants <- PlotIdentifiedModelVariantSignatures(df.sigs.probs = df.curr.variant.signatures.model,
                                                               df.signatures.colors = GetDbsSignaturesColors())
        # Signature by variant plots
        if (length(plots.variants$signature.group.plots) > 0) { # signature group
          for(j in 1:length(plots.variants$signature.group.plots)) {
            ggsave(plot = plots.variants$signature.group.plots[[j]],
                   filename = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                                     sample.id, "_", signatures.count, "_Signature_Model_Variant_Signature_Groups_",
                                     plots.variants$signature.group.plots.names[j], ".svg"),
                   width = 16, height = 9, dpi = 600)
          }
        }
        if (length(plots.variants$signature.plots) > 0) { # signature
          for(j in 1:length(plots.variants$signature.plots)) {
            ggsave(plot = plots.variants$signature.plots[[j]],
                   filename = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                                     sample.id, "_", signatures.count, "_Signature_Model_Variant_Signatures_",
                                     plots.variants$signature.plots.names[j], ".svg"),
                   width = 16, height = 9, dpi = 600)
          }
        }
      }
    }

    # Most parsimonious model
    # Signatures
    plots.models <- PlotIdentifiedModel(identified.model = results[1,],
                                        reference = reference,
                                        plot.theme = plot.theme)
    ggsave(plot = plots.models$plot.merged,
           filename = paste0(save.dir, sample.id, "_Best_Model.svg"),
           width = 16, height = 10, dpi = 600)
    write.table(x = results[1,],
                file = paste0(save.dir, sample.id, "_Best_Model.tsv"),
                sep = "\t", row.names = FALSE)

    plot.weights <- PlotIdentifiedModelSignatureWeights(identified.model = results[1,],
                                                        df.signatures.colors = GetDbsSignaturesColors(),
                                                        plot.theme = plot.theme)
    ggsave(plot = plot.weights,
           filename = paste0(save.dir, sample.id, "_Best_Model_Weights.svg"),
           width = 16, height = 10, dpi = 600)

    # Signatures by variant
    if (analyze.variants == TRUE) {
      best.model.signatures.count <- as.integer(results[1, 'Signatures_Count'])
      best.model.variant.signatures <- results.variants[results.variants$Signatures_Count == best.model.signatures.count,]
      plots.variants <- PlotIdentifiedModelVariantSignatures(df.sigs.probs = best.model.variant.signatures,
                                                             df.signatures.colors = GetDbsSignaturesColors())
      if (length(plots.variants$signature.group.plots) > 0) { # signature group
        for(j in 1:length(plots.variants$signature.group.plots)) {
          ggsave(plot = plots.variants$signature.group.plots[[j]],
                 filename = paste0(save.dir, sample.id, "_Best_Model_Variant_Signature_Groups_",
                                   plots.variants$signature.group.plots.names[j], ".svg"),
                 width = 16, height = 9, dpi = 600)
        }
      }
      if (length(plots.variants$signature.plots) > 0) { # signature
        for(j in 1:length(plots.variants$signature.plots)) {
          ggsave(plot = plots.variants$signature.plots[[j]],
                 filename = paste0(save.dir, sample.id, "_Best_Model_Variant_Signatures_",
                                   plots.variants$signature.plots.names[j], ".svg"),
                 width = 16, height = 9, dpi = 600)
        }
      }

      write.table(x = best.model.variant.signatures,
                  file = paste0(save.dir, sample.id, "_Best_Model_Variant_Signatures.tsv"),
                  sep = "\t", row.names = FALSE)
    }
  } else {
    PrintLog("Step 6/6 - SKIP saving data.")
  }

  return(list(results = results,
              results.variants = results.variants))
}

#' @title Identify indel mutational signatures
#'
#' @description
#' This function identifies indel mutational signatures.
#'
#' @param input Either VCF file path or data.frame. If the input is a data.frame, it must include the following columns: Chr, Pos, Ref, Alt.
#' @param bsg BSgenome object.
#' @param sample.id Sample ID that will be used to name output files (default: 'Sample').
#' @param reference A data.frame with the following columns: Mutation_Type, and names of indel signatures
#' (default: a data.fram returned from \code{\link{GetPcawgIndelSignaturesData(version = "SigProfiler")}}).
#' @param target.signatures Signatures to be considered for identification (default: an array returned from \code{\link{GetPcawgIndelSignaturesNames(version = "SigProfiler")}}).
#' @param plot.theme A data.frame returned from \code{\link{GetIndelSignaturesPlotTheme}}.
#' @param analyze.variants.column.gene Name of column in the data.frame corresponding to the gene name or ID (e.g. "Gene.refGene" if using ANNOVAR data).
#' @param analyze.variants.column.group Name of column in the data.frame corresponding to the variant group for plotting purposes (e.g. "Func.refGene" if using ANNOVAR data).
#' @param analyze.variants A boolean value that indicates whether variant-level signature analysis should be performed (default: FALSE).
#' @param n.cores Number of cores to use.
#' @param max.len Maximum number of bases allowed for a small insertion and deletion (indels bigger than this will be excluded; default: 25).
#' @param padding.len Number of bases to use for upstream and downstream sequences (default: 80).
#' @param combn.m Number of signatures to consider in each step. 'm' parameter in combn function (default: 3).
#' @param n.max.signatures Maximum number of signatures to identify. Recommended: n.max.signatures >= initial.exploration.combn.m (default: 7).
#' @param min.probability Minimum probability to attribute to a signature (default: 0.01).
#' @param zeta.value A float value that is added to the data frequency (default: 1e-10).
#' @param save Save resulting files if TRUE, otherwise do not save (default: TRUE).
#' @param save.dir Save directory path (default: NULL).
#'
#' @return A list with the following elements:
#' @return results: a data.frame with the following columns:
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
#' @return results.variants: a data.frame with the following columns:
#'
#' @export
#' @import dplyr
IdentifyIndelSignatures <- function(input,
                                    bsg,
                                    sample.id = "Sample",
                                    reference = GetPcawgIndelSignaturesData(version = "SigProfiler"),
                                    target.signatures = GetPcawgIndelSignaturesNames(version = "SigProfiler"),
                                    plot.theme = GetIndelSignaturesPlotTheme(),
                                    analyze.variants.column.gene,
                                    analyze.variants.column.group,
                                    analyze.variants = FALSE,
                                    n.cores = 2,
                                    max.len = 25,
                                    padding.len = 80,
                                    combn.m = 3,
                                    n.max.signatures = 7,
                                    min.probability = 0.01,
                                    zeta.value = 1e-10,
                                    save = TRUE,
                                    save.dir = NULL) {
  # Step 1. Check integrity of inputs
  PrintLog("Step 1/6 - Checking integrity of input parameters.")
  if (is.data.frame(input) == TRUE) { # input is data.frame
    input.type <- "data.frame"
  } else {
    if (file.exists(input) == TRUE) { # input is VCF file
      input.type <- "file"
    } else {
      PrintLog("The parameter 'input' appears to be a file but does not exist.", type = "ERROR")
      return()
    }
    if (analyze.variants == TRUE) {
      PrintLog("The parameter 'analyze.variants' cannot be TRUE if the input is a VCF file", type = "ERROR")
      return()
    }
  }
  if (save == TRUE) {
    if (is.null(save.dir) == TRUE) {
      PrintLog("The parameter 'save.dir' must be set if the parameter 'save' is TRUE.", type = "ERROR")
      return()
    } else {
      if (identical(strsplit(x = save.dir, split = "")[[1]], "/") == FALSE) {
        save.dir <- paste0(save.dir, "/")
      }
      if (dir.exists(save.dir) == FALSE) {
        dir.create(path = save.dir, recursive = TRUE, showWarnings = FALSE)
      }
    }
  }
  if (analyze.variants == TRUE) {
    if (is.null(analyze.variants.column.gene)) {
      PrintLog("The parameter 'analyze.variants.column.gene' must not be null if the parameter 'analyze.variants' is TRUE.", type = "ERROR")
      return()
    }
    if (is.null(analyze.variants.column.group)) {
      PrintLog("The parameter 'analyze.variants.column.group' must not be null if the parameter 'analyze.variants' is TRUE.", type = "ERROR")
      return()
    }
  }
  if (is.numeric(n.cores) == FALSE) {
    PrintLog("The parameter 'n.cores' must be a number.", type = "ERROR")
    return()
  } else {
    if (n.cores <= 0) {
      PrintLog("The parameter 'n.cores' must be a positive number", type = "ERROR")
      return()
    }
  }
  if (is.numeric(max.len) == FALSE) {
    PrintLog("The parameter 'max.len' must be a number.", type = "ERROR")
    return()
  }
  if (is.numeric(padding.len) == FALSE) {
    PrintLog("The parameter 'padding.len' must be a number.", type = "ERROR")
    return()
  }
  if (is.numeric(combn.m) == FALSE) {
    PrintLog("The parameter 'combn.m' must be a number.", type = "ERROR")
    return()
  }
  if (is.numeric(n.max.signatures) == FALSE) {
    PrintLog("The parameter 'n.max.signatures' must be a number.", type = "ERROR")
    return()
  }
  if (is.numeric(min.probability) == FALSE) {
    PrintLog("The parameter 'min.probability' must be a number.", type = "ERROR")
    return()
  }
  if (is.numeric(zeta.value) == FALSE) {
    PrintLog("The parameter 'zeta.value' must be a number.", type = "ERROR")
    return()
  }

  # Step 2. Preprocess input data
  PrintLog("Step 2/6 - Preprocessing input variants.")
  if (analyze.variants) {
    input$Variant_Gene <- input[[`analyze.variants.column.gene`]]
    input$Variant_Group <- input[[`analyze.variants.column.group`]]
  }
  if (input.type == "data.frame") {
    prepared.indel.data <- PrepareIndelDataFrame(df = input,
                                                 bsg = bsg,
                                                 max.len = max.len,
                                                 padding.len = padding.len,
                                                 analyze.variants = analyze.variants)
  }
  if (input.type == "file") {
    prepared.indel.data <- PrepareIndelVcfFile(vcf.file = input,
                                               bsg = bsg,
                                               max.len = max.len,
                                               padding.len = padding.len)
  }

  # Step 3. Write prepared indel data to files
  if (save == TRUE) {
    PrintLog("Step 3/6 - Writing preprocessed INDEL data to files.")
    write.table(x = prepared.indel.data$df.indel,
                file = paste0(save.dir, sample.id, "_Indels_Data.tsv"),
                sep = "\t", row.names = FALSE)
    write.table(x = prepared.indel.data$df.indel.frequencies,
                file = paste0(save.dir, sample.id, "_Indels_Frequency.tsv"),
                sep = "\t", row.names = FALSE)

    # Create sub-directory for the i-th signature model
    for (i in 1:n.max.signatures) {
      dir.create(paste0(save.dir, sample.id, "_", i, "_Signature_Model/"))
    }
  } else {
    PrintLog("Step 3/6 - SKIP writing preprocessed INDEL data to files.")
  }

  # Step 4. Identify signatures
  PrintLog("Step 4/6 - Identifying INDEL mutational signatures.")
  results <- IdentifySignatures(data = prepared.indel.data$df.indel.frequencies,
                                reference = reference,
                                target.signatures = target.signatures,
                                n.cores = n.cores,
                                combn.m = combn.m,
                                n.max.signatures = n.max.signatures,
                                min.probability = min.probability,
                                zeta.value = zeta.value)
  results$Sample_ID <- sample.id

  # Step 5. Identify signatures per variant
  results.variants <- data.frame()
  if (analyze.variants == TRUE) {
    PrintLog("Step 5/6 - Identifying INDEL mutational signatures for each variant.")
    for (i in 1:nrow(results)) {
      signatures.count <- as.integer(results[i, 'Signatures_Count'])
      df.curr.model.variant.signatures <- IdentifyVariantSignatures(prepared.data = prepared.indel.data,
                                                                    identified.model = results[i,],
                                                                    reference = reference)
      df.curr.model.variant.signatures$Signatures_Count <- signatures.count
      results.variants <- bind_rows(results.variants, df.curr.model.variant.signatures)
    }
    PrintLog("Finished running indel signature identification by variant.")
  } else {
    PrintLog("Step 5/6 - SKIP identification of SBS mutational signatures for each variant.")
  }

  # Step 6. Save data
  if (save == TRUE) {
    PrintLog("Step 6/6 - Saving data.")
    write.table(x = results,
                file = paste0(save.dir, sample.id, "_All_Models.tsv"),
                sep = "\t", row.names = FALSE)

    if (analyze.variants == TRUE) {
      write.table(x = results.variants,
                  file = paste0(save.dir, sample.id, "_All_Models_Variant_Signatures.tsv"),
                  sep = "\t", row.names = FALSE)
    }

    # All models
    for (i in 1:nrow(results)) {
      # Signatures
      signatures.count <- as.integer(results[i, 'Signatures_Count'])
      curr.model.plots <- PlotIdentifiedModel(identified.model = results[i,],
                                              reference = reference,
                                              plot.theme = plot.theme)
      ggsave(plot = curr.model.plots$plot.merged,
             filename = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                               sample.id, "_", signatures.count, "_Signature_Model.svg"),
             width = 16, height = 10, dpi = 600)

      curr.model.weight.plot <- PlotIdentifiedModelSignatureWeights(identified.model = results[i,],
                                                                    df.signatures.colors = GetIndelSignaturesColors(),
                                                                    plot.theme = plot.theme)
      ggsave(plot = curr.model.weight.plot,
             filename = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                               sample.id, "_", signatures.count, "_Signature_Model_Weights.svg"),
             width = 16, height = 10, dpi = 600)

      # Signatures by variant
      if (analyze.variants == TRUE) {
        df.curr.variant.signatures.model <- results.variants[results.variants$Signatures_Count == signatures.count,]
        write.table(x = df.curr.variant.signatures.model,
                    file = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                                  sample.id, "_", signatures.count, "_Signature_Model_Variant_Signatures.tsv"),
                    sep = "\t", row.names = FALSE)
        plots.variants <- PlotIdentifiedModelVariantSignatures(df.sigs.probs = df.curr.variant.signatures.model,
                                                               df.signatures.colors = GetIndelSignaturesColors())
        # Signature by variant plots
        if (length(plots.variants$signature.group.plots) > 0) { # signature group
          for(j in 1:length(plots.variants$signature.group.plots)) {
            ggsave(plot = plots.variants$signature.group.plots[[j]],
                   filename = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                                     sample.id, "_", signatures.count, "_Signature_Model_Variant_Signature_Groups_",
                                     plots.variants$signature.group.plots.names[j], ".svg"),
                   width = 16, height = 9, dpi = 600)
          }
        }
        if (length(plots.variants$signature.plots) > 0) { # signature
          for(j in 1:length(plots.variants$signature.plots)) {
            ggsave(plot = plots.variants$signature.plots[[j]],
                   filename = paste0(save.dir, sample.id, "_", signatures.count, "_Signature_Model/",
                                     sample.id, "_", signatures.count, "_Signature_Model_Variant_Signatures_",
                                     plots.variants$signature.plots.names[j], ".svg"),
                   width = 16, height = 9, dpi = 600)
          }
        }
      }
    }

    # Most parsimonious model
    # Signatures
    plots.models <- PlotIdentifiedModel(identified.model = results[1,],
                                        reference = reference,
                                        plot.theme = plot.theme)
    ggsave(plot = plots.models$plot.merged,
           filename = paste0(save.dir, sample.id, "_Best_Model.svg"),
           width = 16, height = 10, dpi = 600)
    write.table(x = results[1,],
                file = paste0(save.dir, sample.id, "_Best_Model.tsv"),
                sep = "\t", row.names = FALSE)

    plot.weights <- PlotIdentifiedModelSignatureWeights(identified.model = results[1,],
                                                        df.signatures.colors = GetIndelSignaturesColors(),
                                                        plot.theme = plot.theme)
    ggsave(plot = plot.weights,
           filename = paste0(save.dir, sample.id, "_Best_Model_Weights.svg"),
           width = 16, height = 10, dpi = 600)

    # Signatures by variant
    if (analyze.variants == TRUE) {
      best.model.signatures.count <- as.integer(results[1, 'Signatures_Count'])
      best.model.variant.signatures <- results.variants[results.variants$Signatures_Count == best.model.signatures.count,]
      plots.variants <- PlotIdentifiedModelVariantSignatures(df.sigs.probs = best.model.variant.signatures,
                                                             df.signatures.colors = GetIndelSignaturesColors())
      if (length(plots.variants$signature.group.plots) > 0) { # signature group
        for(j in 1:length(plots.variants$signature.group.plots)) {
          ggsave(plot = plots.variants$signature.group.plots[[j]],
                 filename = paste0(save.dir, sample.id, "_Best_Model_Variant_Signature_Groups_",
                                   plots.variants$signature.group.plots.names[j], ".svg"),
                 width = 16, height = 9, dpi = 600)
        }
      }
      if (length(plots.variants$signature.plots) > 0) { # signature
        for(j in 1:length(plots.variants$signature.plots)) {
          ggsave(plot = plots.variants$signature.plots[[j]],
                 filename = paste0(save.dir, sample.id, "_Best_Model_Variant_Signatures_",
                                   plots.variants$signature.plots.names[j], ".svg"),
                 width = 16, height = 9, dpi = 600)
        }
      }

      write.table(x = best.model.variant.signatures,
                  file = paste0(save.dir, sample.id, "_Best_Model_Variant_Signatures.tsv"),
                  sep = "\t", row.names = FALSE)
    }
  } else {
    PrintLog("Step 6/6 - SKIP saving data.")
  }

  return(list(results = results,
              results.variants = results.variants))
}
