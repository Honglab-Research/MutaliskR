# Plotting functions.
#
# Author: Andy Jinseok Lee


#' @title Plots spectrum
#'
#' @description This function plots spectrum of mutation frequencies.
#'
#' @param y A numeric vector corresponding to the y-axis values.
#' @param plot.theme Plot theme (returned from GetSbsSignaturesPlotTheme).
#'
#' @return A list with the following items:
#' \item{plot}{plot.}
#' \item{legend}{legend plot.}
#'
#' @export
#' @import ggplot2
#' @import ggpubr
#' @importFrom plyr mapvalues
PlotSpectrum <- function(y, plot.theme) {
  names(plot.theme$strip.labels) <- plot.theme$legend.labels

  # Prepare data.frame to plot
  df <- data.frame(x = plot.theme$x.axis.labels,
                   y = y,
                   group = plot.theme$groups,
                   stringsAsFactors = FALSE)
  df$group <- plyr::mapvalues(df$group, from = unique(df$group), to = plot.theme$legend.labels)
  df$group = factor(df$group, levels = plot.theme$legend.labels)

  # Plot
  p <- ggplot(data = df, aes(x = x, y = y, fill = group)) +
    geom_bar(stat = "identity") +
    xlab("") + ylab(plot.theme$y.axis.title) + labs(fill = plot.theme$legend.title) +
    scale_fill_manual(values = plot.theme$strip.colors) +
    scale_y_continuous(expand = c(0, 0), limits = c(plot.theme$y.axis.min, plot.theme$y.axis.max)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = plot.theme$x.axis.text.size,
                                     angle = plot.theme$x.axis.text.angle,
                                     hjust = plot.theme$x.axis.text.hjust,
                                     vjust = plot.theme$x.axis.text.vjust),
          axis.text.y = element_text(size = plot.theme$y.axis.text.size),
          axis.title.y = element_text(size = plot.theme$y.axis.title.size),
          legend.title = element_text(size = plot.theme$legend.title.size),
          legend.text = element_text(size = plot.theme$legend.text.size),
          strip.background = element_rect(fill = "white", colour = "white"),
          strip.text.x = element_text(size = plot.theme$strip.text.size)) +
    facet_grid(. ~ group,
               scales = "free_x",
               space = "free_x",
               labeller = labeller(group = plot.theme$strip.labels))

  # Get legend
  legend <- get_legend(p)

  # Remove legend
  p <- p + theme(legend.position = "none")

  # Assign color to top strip
  g <- ggplot_gtable(ggplot_build(p))
  strip.right <- which(grepl('strip-t', g$layout$name))

  k <- 1
  for (i in strip.right) {
    # Background
    j <- which(grepl('background', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- plot.theme$strip.colors[k]

    # Text
    h <- which(grepl('titleGrob', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[h]]$children[[1]]$gp$col <- plot.theme$strip.text.colors[k]

    k <- k + 1
  }

  plot <- as_ggplot(x = g)
  return(list(plot = plot,
              legend = legend))
}

#' @title Plots signature weights of an identified model
#'
#' @description This function plots signature weights of an identified model.
#'
#' @param identified.model A row of the data.frame returned from the function IdentifySignatures.
#' @param df.signature.colors A data.frame with the following columns: Mutation_Type, and names of signatures.
#'
#' @return ggplot object.
#'
#' @export
#' @import ggplot2
#' @import ggpubr
#' @import stringr
#' @importFrom plyr mapvalues
PlotIdentifiedModelSignatureWeights <- function(identified.model,
                                                df.signatures.colors,
                                                plot.theme) {
  signatures <- as.character(str_split(string = identified.model$Signatures, pattern = "\\,")[[1]])
  signatures.weights <- as.numeric(str_split(string = identified.model$Signatures_Weights, pattern = "\\,")[[1]])

  df.plot <- data.frame(Signature = signatures,
                        Signature_Weight = signatures.weights * 100,
                        stringsAsFactors = FALSE)

  # Signature
  fill.colors <- df.signatures.colors$Color
  names(fill.colors) <- df.signatures.colors$Signature

  df.plot <- df.plot[order(+df.plot$Signature_Weight),]
  df.plot$Signature <- factor(x = df.plot$Signature, levels = unique(df.plot$Signature))
  g <- ggplot(df.plot, aes(x = Signature, y = Signature_Weight, fill = Signature, label = sprintf("%0.1f", round(Signature_Weight, digits = 1)))) +
    geom_bar(position = "stack", stat = "identity") +
    geom_text(size = plot.theme$y.axis.text.size / 3, position = position_dodge(width = 1), hjust = -0.5) +
    xlab("Signature") + ylab("Contribution Probability (%)") +
    scale_fill_manual(values = fill.colors) +
    scale_y_continuous(limits = c(0, max(df.plot$Signature_Weight) * 1.1)) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          axis.text.x = element_text(size = plot.theme$y.axis.text.size, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = plot.theme$y.axis.text.size, hjust = 1, vjust = 0.5),
          axis.title.x = element_text(size = plot.theme$y.axis.title.size),
          axis.title.y = element_text(size = plot.theme$y.axis.title.size)) +
    theme(legend.position = "none") +
    coord_flip()

  # Signature group legend
  df.signatures.colors <- df.signatures.colors[df.signatures.colors$Signature %in% unique(df.plot$Signature),]
  df.plot$Signature_Group <- plyr::mapvalues(
    df.plot$Signature,
    from = df.signatures.colors$Signature,
    to = df.signatures.colors$Signature_Group
  )
  fill.colors <- df.signatures.colors$Color
  names(fill.colors) <- df.signatures.colors$Signature_Group
  df.temp <- data.frame(Signature_Group = unique(df.plot$Signature_Group),
                        Signature_Weight = 1)
  p <- ggplot(df.temp, aes(x = Signature_Group, y = Signature_Weight, fill = Signature_Group)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = fill.colors) +
    labs(fill = "Etiology") +
    theme(legend.text = element_text(size = plot.theme$legend.text.size),
          legend.title = element_text(size = plot.theme$legend.title.size))
  legend.plot <- as_ggplot(get_legend(p = p))

  plot <- ggarrange(plotlist = list(g, legend.plot),
                    ncol = 2, widths = c(2,1), align = "v")
  return(plot)
}

#' @title Plots identified model
#'
#' @description This function plots an identified model.
#'
#' @param identified.model A row of the data.frame returned from the function IdentifySignatures.
#' @param reference A data.frame with the following columns: Mutation_Type, and names of signatures.
#' @param plot.theme A data.frame returned from \code{\link{GetSbsSignaturesPlotTheme}},
#' \code{\link{GetDbsSignaturesPlotTheme}} or \code{\link{GetIndelSignaturesPlotTheme}}.
#'
#' @return A list with the following items:
#' \item{plot.merged}{Merged (plot.observed, plot.reconstructed.spectrum, plot.residual.spectrum, plot.legend) plot.}
#' \item{plot.observed}{Observed mutation frequencies plot.}
#' \item{plot.reconstructed.spectrum}{Reconstructed mutation frequencies plot.}
#' \item{plot.residual.spectrum}{Residual mutation frequencies plot.}
#' \item{plot.legend}{Legend plot.}
#'
#' @export
#' @import ggplot2
#' @import ggpubr
PlotIdentifiedModel <- function(identified.model,
                                reference,
                                plot.theme) {
  # Fetch relevant data
  mutation.types <- as.character(strsplit(as.character(identified.model$Mutation_Types), split = "\\,")[[1]])
  observed.spectrum <- as.numeric(strsplit(as.character(identified.model$Observed_Spectrum), split = "\\,")[[1]])
  reconstructed.spectrum <- as.numeric(strsplit(as.character(identified.model$Reconstructed_Spectrum), split = "\\,")[[1]])
  residual.spectrum <- as.numeric(strsplit(as.character(identified.model$Residual_Spectrum), split = "\\,")[[1]])

  # Determine the y.axis.min.value and y.axis.max.value
  y.axis.min.value <- min(residual.spectrum) * 1.1
  y.axis.max.value <- max(max(observed.spectrum),
                          max(reconstructed.spectrum),
                          max(residual.spectrum))
  y.axis.max.value <- min(1, y.axis.max.value * 1.3)

  # Generate sub plots
  plot.theme$y.axis.min <- 0
  plot.theme$y.axis.max <- y.axis.max.value
  plot.theme$y.axis.title <- "Observed Spectrum Fraction"
  plot.observed.spectrum <- PlotSpectrum(y = observed.spectrum, plot.theme = plot.theme)

  plot.theme$y.axis.min <- 0
  plot.theme$y.axis.max <- y.axis.max.value
  plot.theme$y.axis.title <- "Reconstructed Spectrum Fraction"
  plot.reconstructed.spectrum <- PlotSpectrum(y = reconstructed.spectrum, plot.theme = plot.theme)

  plot.theme$y.axis.min <- y.axis.min.value
  plot.theme$y.axis.max <- y.axis.max.value
  plot.theme$y.axis.title <- "Residual Spectrum Fraction"
  plot.residual.spectrum <- PlotSpectrum(y = residual.spectrum, plot.theme = plot.theme)

  plot.merged <- ggarrange(plotlist = list(plot.observed.spectrum$plot,
                                           plot.reconstructed.spectrum$plot,
                                           plot.residual.spectrum$plot),
                           nrow = 3,
                           legend = "none",
                           align = "hv")

  plot.merged <- ggarrange(plotlist = list(plot.merged,
                                           plot.observed.spectrum$legend),
                           ncol = 2,
                           widths = c(4,1))

  return(list(plot.merged = plot.merged,
              plot.observed = plot.observed.spectrum,
              plot.reconstructed.spectrum = plot.reconstructed.spectrum,
              plot.residual.spectrum = plot.residual.spectrum,
              plot.legend = plot.observed.spectrum$legend))
}
