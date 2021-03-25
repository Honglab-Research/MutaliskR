# Plotting function for a group of samples.
#
# Author: Andy Jinseok Lee


#' @title Sorts a data.frame before plotting a stacked barplot
#'
#' @description This function sorts a data.frame before plotting a stacked barplot.
#'
#' @param df A data.frame.
#' @param ordered.features An ordered character vector of columns to sort in df.
#' @param x.axis.var String value for x axis variable.
#'
#' @return A sorted data.frame.
#'
#' @export
SortDataFrameForStackedBarPlot <- function(df, ordered.features, x.axis.var) {
  sort.string <- c()
  for (curr.feature in ordered.features) {
    sort.string <- c(sort.string, paste0("+df$`", curr.feature, "`"))
  }
  sort.string <- paste0(sort.string, collapse = ",")
  sort.string <- paste0("df[order(", sort.string, "),]")
  df <- eval(parse(text = sort.string))
  return(df)
}

#' @title Plots identified models
#'
#' @description This function plots identified models for a group of samples IDs.
#'
#' @param df.models A data.frame of all best models appended by rbind with Sample_ID added as a column.
#' @param df.signatures.colors A data.frame returned from
#' \code{\link{GetSbsSignaturesColors}}, \code{\link{GetDbsSignaturesColors}}, or \code{\link{GetIndelSignaturesColors}}.
#' @param signature.group.order A character vector specifying the order of the signature groups.
#' @param x.axis.text.size A numeric value that sets the x.axis.text size.
#' @param y.axis.text.size A numeric value that sets the y.axis.text size.
#' @param y.axis.title.size A numeric value that sets the y.axis.title size.
#' @param legend.title.size A numeric value that sets the legend title size.
#' @param legend.text.size A numeric value that sets the legend text size.
#'
#' @return A ggplot object.
#'
#' @export
#' @import ggplot2
#' @import dplyr
#' @importFrom plyr mapvalues
#' @import scales
PlotIdentifiedModels <- function(df.models,
                                 df.signatures.colors,
                                 signature.group.order = c(),
                                 sample.ids.order = c(),
                                 apply.facet = FALSE,
                                 df.sample.ids.groups = c(),
                                 sample.ids.groups.order = c(),
                                 multiply.by.mutations.count = FALSE,
                                 y.axis.max = NULL,
                                 plot.title = "",
                                 title.size = 12,
                                 x.axis.text.size = 8,
                                 y.axis.text.size = 10,
                                 y.axis.title.size = 12,
                                 legend.title.size = 12,
                                 legend.text.size = 10) {
  UnwrapIdentifiedModels <- function(df.models) {
    df.models.unwrapped <- data.frame()
    y.axis.values <- c()
    for (curr.sample.id in unique(df.models$Sample_ID)) {
      df.curr.sample.models <- df.models[df.models$Sample_ID == curr.sample.id,]
      df.curr.sample.models.temp <- data.frame()
      for (i in 1:nrow(df.curr.sample.models)) {
        signatures <- df.curr.sample.models[i,'Signatures']
        signatures <- strsplit(x = signatures, split = "\\,")[[1]]
        signatures.weights <- df.curr.sample.models[i,'Signatures_Weights']
        signatures.weights <- as.numeric(strsplit(x = signatures.weights, split = "\\,")[[1]])
        mutations.count <- as.numeric(df.curr.sample.models[i, "Mutations_Count"])
        if (multiply.by.mutations.count) {
          signatures.weights <- signatures.weights * mutations.count
        } else {
          signatures.weights <- signatures.weights * 100
        }

        df.temp <- data.frame(Signature = signatures,
                              Signature_Weight = signatures.weights,
                              stringsAsFactors = FALSE)
        df.curr.sample.models.temp <- bind_rows(df.curr.sample.models.temp, df.temp)
      }

      # Aggregate
      df.curr.sample.models.aggregated <- df.curr.sample.models.temp %>%
        dplyr::group_by(Signature) %>%
        dplyr::summarise(Signature_Weight = sum(Signature_Weight))
      df.curr.sample.models.aggregated$Sample_ID <- curr.sample.id
      df.models.unwrapped <- bind_rows(df.models.unwrapped, df.curr.sample.models.aggregated)
      y.axis.values <- c(y.axis.values, sum(df.curr.sample.models.aggregated$Signature_Weight))
    }

    return(list(df.models.unwrapped = df.models.unwrapped, y.axis.max = max(y.axis.values)))
  }
  unwrapped.data <- UnwrapIdentifiedModels(df.models = df.models)
  df.models.unwrapped <- unwrapped.data$df.models.unwrapped

  if (is.null(y.axis.max)) {
    y.axis.max <- unwrapped.data$y.axis.max
  }

  # Map signatures to signature groups
  df.models.unwrapped$Signature_Group <- plyr::mapvalues(
    df.models.unwrapped$Signature,
    from = df.signatures.colors$Signature,
    to = df.signatures.colors$Signature_Group
  )

  # Get probability sum
  if (length(signature.group.order) == 0) {
    all.signature.groups <- unique(df.models.unwrapped$Signature_Group)
    df.models.unwrapped.summed <- df.models.unwrapped %>%
      group_by(Signature_Group) %>%
      summarise_each(funs(sum), Signature_Weight)
    df.models.unwrapped.summed <- df.models.unwrapped.summed[order(+df.models.unwrapped.summed$Signature_Weight),]

    # Get signature group order
    signature.group.order <- unique(df.models.unwrapped.summed$Signature_Group)
    df.models.unwrapped$Signature_Group <- factor(x = df.models.unwrapped$Signature_Group,
                                                  levels = signature.group.order)
  } else {
    df.models.unwrapped$Signature_Group <- factor(x = df.models.unwrapped$Signature_Group,
                                                  levels = signature.group.order)
  }

  # Fetch signature colors
  fill.colors <- df.signatures.colors$Color
  names(fill.colors) <- df.signatures.colors$Signature_Group

  # Sort sample IDs
  if (length(sample.ids.order) == 0) {
    df.sorting.info <- data.frame()
    for (curr.sample.id in unique(df.models.unwrapped$Sample_ID)) {
      df.curr.sample.id <- data.frame(Sample_ID = curr.sample.id)

      if (multiply.by.mutations.count) {
        signature.weights.sum <- sum(df.models.unwrapped[df.models.unwrapped$Sample_ID == curr.sample.id,"Signature_Weight"])
        df.temp <- data.frame(Signature_Weights_Sum = signature.weights.sum)
        df.curr.sample.id <- bind_cols(df.curr.sample.id, df.temp)
      }

      for (curr.signature.group in signature.group.order) {
        conditions <- ((df.models.unwrapped$Sample_ID == curr.sample.id) &
                         (df.models.unwrapped$Signature_Group == curr.signature.group))
        df.temp <- df.models.unwrapped[conditions,]
        if (nrow(df.temp) == 0) {
          df.curr.signature.group.weights <- data.frame(a = 0)
        } else {
          df.curr.signature.group.weights <- data.frame(a = sum(df.temp$Signature_Weight))
        }
        colnames(df.curr.signature.group.weights) <- curr.signature.group
        df.curr.sample.id <- bind_cols(df.curr.sample.id, df.curr.signature.group.weights)
      }
      df.sorting.info <- bind_rows(df.sorting.info, df.curr.sample.id)
    }

    ordered.features <- colnames(df.sorting.info)
    ordered.features <- ordered.features[ordered.features != "Sample_ID"]
    df.sorting.info <- SortDataFrameForStackedBarPlot(df = df.sorting.info,
                                                      ordered.features = ordered.features,
                                                      x.axis.var = "Sample_ID")
    df.models.unwrapped$Sample_ID <- factor(x = df.models.unwrapped$Sample_ID,
                                            levels = rev(df.sorting.info$Sample_ID))
  } else {
    df.models.unwrapped$Sample_ID <- factor(x = df.models.unwrapped$Sample_ID,
                                            levels = sample.ids.order)
  }

  # Plot
  if (multiply.by.mutations.count) {
    y.axis.title <- "Signature Mutations Count"
  } else {
    y.axis.title <- "Signature Contribution Probability"
  }

  if (apply.facet) {
    df.models.unwrapped$Group <- plyr::mapvalues(
      x = df.models.unwrapped$Sample_ID,
      from = df.sample.ids.groups$Sample_ID,
      to = df.sample.ids.groups$Group
    )
    df.models.unwrapped$Group <- factor(x = df.models.unwrapped$Group,
                                        levels = sample.ids.groups.order)
    g <- ggplot(data = df.models.unwrapped,
                aes(x = Sample_ID,
                    y = Signature_Weight,
                    fill = Signature_Group)) +
      geom_bar(position = "stack", stat = "identity") +
      ggtitle(plot.title) +
      xlab("Signature Weight (%)") + ylab(y.axis.title) + labs(fill = "Etiology") +
      scale_y_continuous(expand = c(0,0), limits = c(0, y.axis.max), labels = scales::comma) +
      scale_fill_manual(values = fill.colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(size = x.axis.text.size, angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = y.axis.text.size),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = y.axis.title.size),
            plot.title = element_text(size = title.size, hjust = 0.5),
            panel.grid.major.x = element_blank(),
            legend.title = element_text(size = legend.title.size),
            legend.text = element_text(size = legend.text.size)) +
      facet_grid(.~Group, scales = "free", space = 'free_x')

  } else {
    g <- ggplot(data = df.models.unwrapped,
                aes(x = Sample_ID,
                    y = Signature_Weight,
                    fill = Signature_Group)) +
      geom_bar(position = "stack", stat = "identity") +
      ggtitle(plot.title) +
      xlab("Signature Weight (%)") + ylab(y.axis.title) + labs(fill = "Etiology") +
      scale_y_continuous(expand = c(0,0), limits = c(0, y.axis.max), labels = scales::comma) +
      scale_fill_manual(values = fill.colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(size = x.axis.text.size, angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = y.axis.text.size),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = y.axis.title.size),
            plot.title = element_text(size = title.size, hjust = 0.5),
            panel.grid.major.x = element_blank(),
            legend.title = element_text(size = legend.title.size),
            legend.text = element_text(size = legend.text.size))
  }

  return(list(plot = g, df.plot = df.models.unwrapped))
}

#' @title Plots a stacked barplot of signatures probabilities of variants
#'
#' @description This function plots a stacked barplot of signature probabilities of variants.
#'
#' @return A list with the following items:
#' \item{plots.merged}{Merged rainfall plots.}
#' \item{plots}{Rainfall plot for each chromosome.}
#'
#' @export
#' @import ggplot2
#' @import ggpubr
#' @import scales
#' @import tools
#' @import stringr
#' @importFrom plyr mapvalues
PlotIdentifiedModelVariantSignatures <- function(df.sigs.probs,
                                                 df.signatures.colors) {
  # Prepare plotting data
  df.sigs.probs <- df.sigs.probs[df.sigs.probs$Variant_Gene != "",]
  df.sigs.probs$Probability <- df.sigs.probs$Probability * 100
  df.signatures.colors <- df.signatures.colors[df.signatures.colors$Signature %in% unique(df.sigs.probs$Signature),]

  # Define plot helper function
  PlotHelper <- function(df.plot, fill.colors, title, legend.title) {
    g <- ggplot(data = df.plot, aes(x = Variant_Gene, y = Probability, fill = fill)) +
      geom_bar(position = "stack", stat = "identity") + ggtitle(title) +
      xlab("") + ylab("Cumulative Probability (%)") + labs(fill = legend.title) +
      scale_fill_manual(values = fill.colors) +
      scale_y_continuous(expand = c(0,0), label = comma) +
      theme_minimal() + coord_flip() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            plot.title = element_text(size = 12, vjust = 0.5, hjust = 0.5),
            axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 0.5),
            axis.text.y = element_text(size = 8, vjust = 0.5, hjust = 1),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_blank(),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10))
    return(g)
  }

  # Plot type 1. Signature group
  PrepareSignatureGroupPlotData <- function(df.sigs.probs.selected) {
    # Get probability sum
    all.signature.groups <- unique(df.sigs.probs.selected$Signature_Group)
    df.sigs.probs.selected.summed <- df.sigs.probs.selected %>%
      group_by(Signature_Group) %>%
      summarise_each(funs(sum), Probability)
    df.sigs.probs.selected.summed <- df.sigs.probs.selected.summed[order(-df.sigs.probs.selected.summed$Probability),]

    # Get signature group ranking
    signature.group.ranking <- unique(df.sigs.probs.selected.summed$Signature_Group)
    df.sigs.probs.selected$Signature_Group <- factor(x = df.sigs.probs.selected$Signature_Group, levels = signature.group.ranking)

    # Sort by cumulative probability and top signature probability
    AggregateSortInfo <- function(x) {
      curr.gene.name <- as.character(x[['Variant_Gene']])
      df.matched <- df.sigs.probs.selected %>%
        filter(Variant_Gene == curr.gene.name)
      probability.sum <- sum(df.matched$Probability)
      df.return <- data.frame(Variant_Gene = curr.gene.name,
                              Probability_Sum = round(x = probability.sum, digits = 0),
                              stringsAsFactors = FALSE)
      for (curr.signature.group in signature.group.ranking) {
        curr.signature.group.prob <- sum(df.matched[df.matched$Signature_Group == curr.signature.group, "Probability"])
        df.curr.signature.group.prob <- data.frame(a = curr.signature.group.prob)
        colnames(df.curr.signature.group.prob) <- curr.signature.group
        df.return <- bind_cols(df.return, df.curr.signature.group.prob)
      }
      return(df.return)
    }
    df.sorting.info <- apply(X = data.frame(Variant_Gene = unique(df.sigs.probs.selected$Variant_Gene)),
                             MARGIN = 1, FUN = AggregateSortInfo)
    df.sorting.info <- bind_rows(df.sorting.info)
    ordered.features <- colnames(df.sorting.info)
    ordered.features <- ordered.features[ordered.features != "Variant_Gene"]
    df.sorting.info <- SortDataFrameForStackedBarPlot(df = df.sorting.info,
                                                      ordered.features = ordered.features,
                                                      x.axis.var = "Variant_Gene")
    df.sigs.probs.selected$Variant_Gene <- factor(x = df.sigs.probs.selected$Variant_Gene,
                                                  levels = rev(df.sorting.info$Variant_Gene))
    return(df.sigs.probs.selected)
  }

  fill.colors <- df.signatures.colors$Color
  names(fill.colors) <- df.signatures.colors$Signature_Group
  df.sigs.probs$Signature_Group <- mapvalues(df.sigs.probs$Signature,
                                             from = df.signatures.colors$Signature,
                                             to = df.signatures.colors$Signature_Group)

  signature.group.plots <- list()
  signature.group.plots.names <- c()
  plots.idx <- 1
  for (curr.variant.group in unique(df.sigs.probs$Variant_Group)) {
    df.sigs.probs.selected <- df.sigs.probs[df.sigs.probs$Variant_Group == curr.variant.group,]
    if (curr.variant.group == "") { # intergenic variants
      next
    } else if (curr.variant.group == "unknown") {
      next
    } else {
      df.sigs.probs.selected <- PrepareSignatureGroupPlotData(df.sigs.probs.selected = df.sigs.probs.selected)
      df.sigs.probs.selected$fill <- df.sigs.probs.selected$Signature_Group
      title <- paste0(tools::toTitleCase(curr.variant.group), " Variants")
      signature.group.plots[[plots.idx]] <- PlotHelper(df.plot = df.sigs.probs.selected,
                                                       fill.colors = fill.colors,
                                                       title = title,
                                                       legend.title = "Signature Group")
      signature.group.plots.names <- c(signature.group.plots.names, str_replace(string = title, pattern = " ", replacement = "_"))
      plots.idx <- plots.idx + 1
    }
  }

  # Plot type 2. Signature
  PrepareSignaturePlotData <- function(df.sigs.probs.selected) {
    # Get probability sum
    all.signatures <- unique(df.sigs.probs.selected$Signature)
    df.sigs.probs.selected.summed <- df.sigs.probs.selected %>%
      group_by(Signature) %>%
      summarise_each(funs(sum), Probability)
    df.sigs.probs.selected.summed <- df.sigs.probs.selected.summed[order(-df.sigs.probs.selected.summed$Probability),]

    # Get signature ranking
    signature.ranking <- unique(df.sigs.probs.selected.summed$Signature)
    df.sigs.probs.selected$Signature <- factor(x = df.sigs.probs.selected$Signature, levels = signature.ranking)

    # Sort by cumulative probability and top signature probability
    AggregateSortInfo <- function(x) {
      curr.gene.name <- as.character(x[['Variant_Gene']])
      df.matched <- df.sigs.probs.selected %>%
        filter(Variant_Gene == curr.gene.name)
      probability.sum <- sum(df.matched$Probability)
      df.return <- data.frame(Variant_Gene = curr.gene.name,
                              Probability_Sum = round(x = probability.sum, digits = 0),
                              stringsAsFactors = FALSE)
      for (curr.signature in signature.ranking) {
        curr.signature.prob <- sum(df.matched[df.matched$Signature == curr.signature, "Probability"])
        df.curr.signature.prob <- data.frame(a = curr.signature.prob)
        colnames(df.curr.signature.prob) <- curr.signature
        df.return <- bind_cols(df.return, df.curr.signature.prob)
      }
      return(df.return)
    }
    df.sorting.info <- apply(X = data.frame(Variant_Gene = unique(df.sigs.probs.selected$Variant_Gene)),
                             MARGIN = 1, FUN = AggregateSortInfo)
    df.sorting.info <- bind_rows(df.sorting.info)
    ordered.features <- colnames(df.sorting.info)
    ordered.features <- ordered.features[ordered.features != "Variant_Gene"]
    df.sorting.info <- SortDataFrameForStackedBarPlot(df = df.sorting.info,
                                                      ordered.features = ordered.features,
                                                      x.axis.var = "Variant_Gene")
    df.sigs.probs.selected$Variant_Gene <- factor(x = df.sigs.probs.selected$Variant_Gene,
                                                  levels = rev(df.sorting.info$Variant_Gene))
    return(df.sigs.probs.selected)
  }

  fill.colors <- df.signatures.colors$Color
  names(fill.colors) <- df.signatures.colors$Signature

  signature.plots <- list()
  signature.plots.names <- c()
  plots.idx <- 1
  for (curr.variant.group in unique(df.sigs.probs$Variant_Group)) {
    df.sigs.probs.selected <- df.sigs.probs[df.sigs.probs$Variant_Group == curr.variant.group,]
    if (curr.variant.group == "") { # intergenic variants
      next
    } else {
      df.sigs.probs.selected <- PrepareSignaturePlotData(df.sigs.probs.selected = df.sigs.probs.selected)
      df.sigs.probs.selected$fill <- df.sigs.probs.selected$Signature
      title <- paste0(tools::toTitleCase(curr.variant.group), " Variants")
      signature.plots[[plots.idx]] <- PlotHelper(df.plot = df.sigs.probs.selected,
                                                 fill.colors = fill.colors,
                                                 title = title,
                                                 legend.title = "Signature")
      signature.plots.names <- c(signature.plots.names, str_replace(string = title, pattern = " ", replacement = "_"))
      plots.idx <- plots.idx + 1
    }
  }

  return(list(signature.group.plots = signature.group.plots,
              signature.group.plots.names = signature.group.plots.names,
              signature.plots = signature.plots,
              signature.plots.names = signature.plots.names))
}

