# Plotting themes.
#
# Author: Andy Jinseok Lee


#' @title Fetches SBS plotting theme
#'
#' @description This function returns the SBS plotting theme elements.
#'
#' @param x.axis.labels A character vector used to label the x-axis.
#' @param x.axis.text.size x-axis text size.
#' @param x.axis.text.angle x-axis text angle.
#' @param y.axis.text.size y-axis text size.
#' @param y.axis.title.size y-axis title size.
#' @param y.axis.min y-axis lower-bound (minimum range for limits).
#' @param y.axis.max y-axis upper-bound (maximum range for limits).
#' @param y.axis.title y-axis title text.
#' @param groups A character vector of group variable values.
#' @param strip.labels A character vector used to label the strip.
#' @param strip.colors A character vector of hex values used to color the strip.
#' @param strip.text.colors A character vector of hex values used to color the strip text.
#' @param strip.text.size Strip (facet) text size.
#' @param legend.title Legend title text.
#' @param legend.title.size Legend title text size.
#' @param legend.labels A character vector used to label the legend items.
#'
#' @return A data.frame of the SBS signature plot theme elements.
#'
#' @export
#' @import ggplot2
#' @import ggpubr
GetSbsSignaturesPlotTheme <- function(x.axis.labels = PLOT.SBS.X.AXIS.LABELS,
                                      x.axis.text.size = 10,
                                      x.axis.text.angle = 90,
                                      x.axis.text.hjust = 0.5,
                                      x.axis.text.vjust = 0.5,
                                      y.axis.text.size = 12,
                                      y.axis.title.size = 12,
                                      y.axis.min = 0,
                                      y.axis.max = 1,
                                      y.axis.title = "Mutation Fraction",
                                      groups = SBS.MUTATION.TYPES.GROUPS,
                                      strip.labels = PLOT.SBS.STRIP.LABELS,
                                      strip.colors = PLOT.SBS.STRIP.COLORS,
                                      strip.text.colors = PLOT.SBS.STRIP.TEXT.COLORS,
                                      strip.text.size = 12,
                                      legend.title = "Mutation Type",
                                      legend.title.size = 12,
                                      legend.text.size = 10,
                                      legend.labels = PLOT.SBS.LEGEND.LABELS) {
  return(list(x.axis.labels = x.axis.labels,
              x.axis.text.size = x.axis.text.size,
              x.axis.text.angle = x.axis.text.angle,
              x.axis.text.hjust = x.axis.text.hjust,
              x.axis.text.vjust = x.axis.text.vjust,
              y.axis.text.size = y.axis.text.size,
              y.axis.title.size = y.axis.title.size,
              y.axis.min = y.axis.min,
              y.axis.max = y.axis.max,
              y.axis.title = y.axis.title,
              groups = groups,
              strip.labels = strip.labels,
              strip.colors = strip.colors,
              strip.text.colors = strip.text.colors,
              strip.text.size = strip.text.size,
              legend.title = legend.title,
              legend.title.size = legend.title.size,
              legend.text.size = legend.text.size,
              legend.labels = legend.labels))
}

#' @title Fetches DBS plotting theme
#'
#' @description This function returns the DBS plotting theme elements.
#'
#' @param x.axis.labels A character vector used to label the x-axis.
#' @param x.axis.text.size x-axis text size.
#' @param x.axis.text.angle x-axis text angle.
#' @param y.axis.text.size y-axis text size.
#' @param y.axis.title.size y-axis title size.
#' @param y.axis.min y-axis lower-bound (minimum range for limits).
#' @param y.axis.max y-axis upper-bound (maximum range for limits).
#' @param y.axis.title y-axis title text.
#' @param groups A character vector of group variable values.
#' @param strip.labels A character vector used to label the strip.
#' @param strip.colors A character vector of hex values used to color the strip.
#' @param strip.text.colors A character vector of hex values used to color the strip text.
#' @param strip.text.size Strip (facet) text size.
#' @param legend.title Legend title text.
#' @param legend.title.size Legend title text size.
#' @param legend.labels A character vector used to label the legend items.
#'
#' @return A data.frame of the DBS signature plot theme elements.
#'
#' @export
#' @import ggplot2
#' @import ggpubr
GetDbsSignaturesPlotTheme <- function(x.axis.labels = PLOT.DBS.X.AXIS.LABELS,
                                      x.axis.text.size = 10,
                                      x.axis.text.angle = 90,
                                      x.axis.text.hjust = 1,
                                      x.axis.text.vjust = 0.5,
                                      y.axis.text.size = 12,
                                      y.axis.title.size = 12,
                                      y.axis.min = 0,
                                      y.axis.max = 1,
                                      y.axis.title = "Mutation Fraction",
                                      groups = DBS.MUTATION.TYPES.GROUPS,
                                      strip.labels = PLOT.DBS.STRIP.LABELS,
                                      strip.colors = PLOT.DBS.STRIP.COLORS,
                                      strip.text.colors = PLOT.DBS.STRIP.TEXT.COLORS,
                                      strip.text.size = 12,
                                      legend.title = "Mutation Type",
                                      legend.title.size = 12,
                                      legend.text.size = 10,
                                      legend.labels = PLOT.DBS.LEGEND.LABELS) {
  return(list(x.axis.labels = x.axis.labels,
              x.axis.text.size = x.axis.text.size,
              x.axis.text.angle = x.axis.text.angle,
              x.axis.text.hjust = x.axis.text.hjust,
              x.axis.text.vjust = x.axis.text.vjust,
              y.axis.text.size = y.axis.text.size,
              y.axis.title.size = y.axis.title.size,
              y.axis.min = y.axis.min,
              y.axis.max = y.axis.max,
              y.axis.title = y.axis.title,
              groups = groups,
              strip.labels = strip.labels,
              strip.colors = strip.colors,
              strip.text.colors = strip.text.colors,
              strip.text.size = strip.text.size,
              legend.title = legend.title,
              legend.title.size = legend.title.size,
              legend.text.size = legend.text.size,
              legend.labels = legend.labels))
}

#' @title Fetches INDEL plotting theme
#'
#' @description This function returns the INDEL plotting theme elements.
#'
#' @param x.axis.labels A character vector used to label the x-axis.
#' @param x.axis.text.size x-axis text size.
#' @param x.axis.text.angle x-axis text angle.
#' @param y.axis.text.size y-axis text size.
#' @param y.axis.title.size y-axis title size.
#' @param y.axis.min y-axis lower-bound (minimum range for limits).
#' @param y.axis.max y-axis upper-bound (maximum range for limits).
#' @param y.axis.title y-axis title text.
#' @param groups A character vector of group variable values.
#' @param strip.labels A character vector used to label the strip.
#' @param strip.colors A character vector of hex values used to color the strip.
#' @param strip.text.colors A character vector of hex values used to color the strip text.
#' @param strip.text.size Strip (facet) text size.
#' @param legend.title Legend title text.
#' @param legend.title.size Legend title text size.
#' @param legend.labels A character vector used to label the legend items.
#'
#' @return A data.frame of the INDEL signature plot theme elements.
#'
#' @export
#' @import ggplot2
#' @import ggpubr
GetIndelSignaturesPlotTheme <- function(x.axis.labels = PLOT.INDEL.X.AXIS.LABELS,
                                        x.axis.text.size = 10,
                                        x.axis.text.angle = 0,
                                        y.axis.text.size = 12,
                                        x.axis.text.hjust = 0.5,
                                        x.axis.text.vjust = 0.5,
                                        y.axis.title.size = 12,
                                        y.axis.min = 0,
                                        y.axis.max = 1,
                                        y.axis.title = "Mutation Fraction",
                                        groups = INDEL.MUTATION.TYPES.GROUPS,
                                        strip.labels = PLOT.INDEL.STRIP.LABELS,
                                        strip.colors = PLOT.INDEL.STRIP.COLORS,
                                        strip.text.colors = PLOT.INDEL.STRIP.TEXT.COLORS,
                                        strip.text.size = 12,
                                        legend.title = "Mutation Type",
                                        legend.title.size = 12,
                                        legend.text.size = 10,
                                        legend.labels = PLOT.INDEL.LEGEND.LABELS) {
  return(list(x.axis.labels = x.axis.labels,
              x.axis.text.size = x.axis.text.size,
              x.axis.text.angle = x.axis.text.angle,
              x.axis.text.hjust = x.axis.text.hjust,
              x.axis.text.vjust = x.axis.text.vjust,
              y.axis.text.size = y.axis.text.size,
              y.axis.title.size = y.axis.title.size,
              y.axis.min = y.axis.min,
              y.axis.max = y.axis.max,
              y.axis.title = y.axis.title,
              groups = groups,
              strip.labels = strip.labels,
              strip.colors = strip.colors,
              strip.text.colors = strip.text.colors,
              strip.text.size = strip.text.size,
              legend.title = legend.title,
              legend.title.size = legend.title.size,
              legend.text.size = legend.text.size,
              legend.labels = legend.labels))
}


