#' Generate a distinct palette for coloring different clusters.
#'
#' Create a distinct palette for coloring different heatmap clusters. The
#' function returns colors for input into \code{ComplexHeatmap:Heatmap()}.
#'
#' @param n an integer describing the number of colors to generate. Required.
#' @param hues a vector of character strings indicating the R colors available
#' from the \code{colors()} function. These will be used as the base colors for
#' the clustering scheme. Different saturations and values (i.e. darkness)
#' will be generated for each hue. Default is \code{c("red", "cyan", "orange",
#' "blue", "yellow", "purple", "green", "magenta")}
#' @param saturation.range a numeric vector of length 2 with values between 0
#' and 1 giving the range of saturation. The default is \code{c(0.25, 1)}.
#' @param value.range a numeric vector of length 2 with values between 0 and 1
#' giving the range of values. The default is \code{c(0.5, 1)}.
#'
#' @return A vector of distinct colors that have been converted to HEX from
#' HSV.
#'
#' @export
#'
#' @examples
#'
#' distinctColors(10)
#'
distinctColors <- function(n, hues = c("red", "cyan", "orange", "blue",
                                       "yellow", "purple", "green", "magenta"),
                           saturation.range = c(0.7, 1),
                           value.range = c(0.7, 1)) {
  
  if (!(all(hues %in% grDevices::colors()))) {
    stop("Only color names listed in the 'color'",
         " function can be used in 'hues'")
  }
  
  ## Convert R colors to RGB and then to HSV color format
  hues.hsv <- grDevices::rgb2hsv(grDevices::col2rgb(hues))
  
  ## Calculate all combination of saturation/value pairs.
  ## Note that low saturation with low value (i.e. high darkness) is too dark
  ## for all hues. Likewise, high saturation with high value (i.e. low darkness)
  ## is hard to distinguish. Therefore, saturation and value are set to be
  ## anticorrelated.
  num.vs <- ceiling(n / length(hues))
  s <- seq(from = saturation.range[1], to = saturation.range[2],
           length = num.vs)
  v <- seq(from = value.range[2], to = value.range[1], length = num.vs)
  
  ## Create all combinations of hues with saturation/value pairs
  new.hsv <- c()
  for (i in seq_len(num.vs)) {
    temp <- rbind(hues.hsv[1, ], s[i], v[i])
    new.hsv <- cbind(new.hsv, temp)
  }
  
  ## Convert to HEX
  col <- grDevices::hsv(new.hsv[1, ], new.hsv[2, ], new.hsv[3, ])
  
  return(col[seq_len(n)])
}

#' Plot a ComplexHeatmap.
#'
#' This function takes an arbitrary dataset as an input
#' and returns a \code{ComplexHeatmap} plot. The function takes arguments listed here as well
#' as any others to be passed on to \code{ComplexHeatmap::Heatmap()}.
#'
#' If both \code{annotationData = NULL} and \code{annotationColNames = NULL},
#' no annotation bar will be drawn on the heatmap.
#' 
#' Code was adapted from the \code{TBSignatureprofiler} R package.
#'
#' @param inputData an input data object. It should either be of the class
#' \code{SummarizedExperiment} and contain the data and
#' annotation data as columns in the colData, or alternatively be of the classes
#' \code{data.frame} or \code{matrix} and contain only the plotting data.
#' Required.
#' @param annotationData a \code{data.frame} or \code{matrix} of annotation
#' data, with one column. Only required if \code{inputData} is a
#' \code{data.frame} or \code{matrix} of plotting data.
#' The row names must equal those of the \code{inputData} column names.
#' Default is \code{NULL}.
#' @param plot_title a character string with the plot title of the heatmap. The
#' default is \code{NULL}.
#' @param name a character string with the name of the data to be displayed.
#' Default is \code{"Input data"}.
#' @param plottingColNames a vector of the column names in \code{colData} that
#' contain the plotting data. Only required if \code{inputData} is a
#' SummarizedExperiment object.
#' @param annotationColNames a vector of the column names in \code{colData} that
#' contain the annotation data. Only required if \code{inputData} is a
#' \code{SummarizedExperiment}. Default is \code{NULL}.
#' @param colList a named \code{list} of named vectors specifying custom color
#' information to
#' pass to \code{ComplexHeatmap::Heatmap()}. The list should have as many
#' elements as there are annotation columns, and each element name should
#' correspond exactly with the name of each annotation column.
#' The colors in the vector elements should be named according to the
#' levels of the factor in that column's annotation data if the annotation
#' is discrete, or it should be produced with \code{circlize::colorRamp2}
#' if the annotation is continuous.
#' By default, \code{ColorBrewer} color sets will be used.
#' See the the parameter \code{colorSets} for additional details.
#' @param scale logical. Setting \code{scale = TRUE} scales the plotting data.
#' The default is \code{FALSE}.
#' @param showColumnNames logical. Setting \code{showColumnNames = TRUE} will
#' show the column names (i.e. sample names) on the heatmap. The default is
#' \code{TRUE}.
#' @param showRowNames logical. Setting \code{showColumnNames = TRUE} will
#' show the row names (i.e. plotting names) on the heatmap. The default is
#' \code{TRUE}.
#' @param colorSets a vector of names listing the color sets in the order
#' that they should be used in creating the heatmap. By default, this function
#' will use the color sets in the order listed in \code{Usage} for annotation
#' information. You may replace the default with the same collection of sets
#' in order that you want to use them, or provide custom color sets with the
#' \code{colList} parameter.
#' @param choose_color a vector of color names to be interpolated for the
#' heatmap gradient, or a \code{colorRamp} function produced by
#' \code{circlize::colorRamp2}. The default is \code{c("blue", "gray95", "red")}.
#' @param split_heatmap a character string either giving the column title of
#' \code{annotationplotting} containing annotation data for which to split
#' the heatmap rows, or \code{"none"} if no split is desired.
#' @param annotationplotting a \code{data.frame} or \code{matrix} with information
#' to be used
#' in splitting the heatmap. The first column should plotting names. The
#' column of annotation information should be specified in \code{split_heatmap.}
#' Other columns will be ignored. The default is \code{sigAnnotData}.
#' @param column_order a vector of character strings indicating the order in
#' which to manually arrange the heatmap columns. Default is \code{NULL},
#' such that column order is automatically determined via clustering.
#' @param ... Additional arguments to be passed to
#' \code{ComplexHeatmap::Heatmap()}.
#'
#' @return A ComplexHeatmap plot.
#' 
#' @author David Jenkins, Aubrey Odom
#' 
#'
#' @export
#'@examples
#' library(SummarizedExperiment)
#' # Generate some artificial data that shows a difference in Zak_RISK_16
#' mat_testdata <- rbind(matrix(c(rnorm(80), rnorm(80) + 5), 16, 10,
#'                              dimnames = list(paste0("Taxon", seq_len(16)),
#'                                              paste0("sample", seq_len(10)))),
#'                       matrix(rnorm(1000), 100, 10,
#'                              dimnames = list(paste0("Taxon0", seq_len(100)),
#'                                              paste0("sample", seq_len(10)))))
#' cov_mat <- data.frame(sample = c(rep("down", 5), rep("up", 5))) |>
#'   magrittr::set_rownames(paste0("sample", seq_len(10)))
#' 
#' # Example using custom colors for the annotation information
#' color2 <- stats::setNames(c("purple", "black"), c("down", "up"))
#' color.list <- list("sample" = color2)
#' 
#' plot_heatmap(
#'   inputData = mat_testdata,
#'   annotationData = cov_mat,
#'   name = "Data",
#'   plot_title = "Example",
#'   plottingColNames = NULL,
#'   annotationColNames = NULL,
#'   colList = color.list,
#'   scale = FALSE,
#'   showColumnNames = TRUE,
#'   showRowNames = FALSE,
#'   colorSets = c("Set1", "Set2", "Set3", "Pastel1", "Pastel2", "Accent", "Dark2",
#'                 "Paired"),
#'   choose_color = c("blue", "gray95", "red"),
#'   split_heatmap = "none",
#'   column_order = NULL
#' )
#'

plot_heatmap <- function(inputData, annotationData = NULL, plot_title = NULL,
                         name = "Input data",
                         plottingColNames,
                         annotationColNames = NULL,
                         colList = list(), scale = FALSE,
                         showColumnNames = TRUE,
                         showRowNames = TRUE, colorSets = c("Set1", "Set2",
                                                            "Set3", "Pastel1", "Pastel2", "Accent", "Dark2",
                                                            "Paired"),
                         choose_color = c("blue", "gray95", "red"),
                         split_heatmap = "none",
                         annotationplotting = NULL,
                         column_order = NULL,
                         ...) {
  if (methods::is(inputData, "SummarizedExperiment")) {
    if (any(duplicated(plottingColNames))) {
      stop("Duplicate plotting column name is not supported.")
    }
    if (!all(plottingColNames %in% colnames(SummarizedExperiment::colData(inputData)))) {
      stop("plotting column name not found in inputData.")
    }
    if (!is.null(annotationColNames)) {
      if (!all(annotationColNames %in% colnames(SummarizedExperiment::colData(inputData)))) {
        stop("Annotation column name not found in inputData.")
      }
      annotationData <- SummarizedExperiment::colData(inputData)[, annotationColNames, drop = FALSE]
      inputData <- SummarizedExperiment::colData(inputData)[, plottingColNames, drop = FALSE]
    }
  } else {
    if (is.null(annotationData)) {
      stop("annotationData must be provided for a data.frame input object.")
    } else if (!is.null(annotationData)) {
      annotationColNames <- colnames(annotationData)
    }
  }
  if (!is.null(annotationData)) {
    if (nrow(annotationData) == nrow(inputData)) {
      if (!all(rownames(annotationData) == rownames(inputData))) {
        stop("Annotation data and plotting data does not match.")
      }
    } else if (nrow(annotationData) == ncol(inputData)) {
      if (!all(rownames(annotationData) == colnames(inputData))) {
        stop("Annotation data and plotting data does not match.")
      }
      inputData <- t(inputData)
    } else {
      stop("Annotation data and plotting data does not match.")
    }
  }
  sigresults <- t(as.matrix(inputData))
  if (scale) {
    sigresults <- t(scale(t(sigresults)))
  }
  # To split heatmap by plotting annotation
  if (split_heatmap != "none") {
    if (!(split_heatmap %in% colnames(annotationplotting))) {
      stop("The column specified in 'split_heatmap' must be in the matrix or data.frame
           provided by 'annotationplotting'")
    }
  }
  ann_data <- annotationplotting[annotationplotting[, 1] %in%
                                    plottingColNames, ]
  if (split_heatmap == "none") {
    row_split_pass <- c()
  } else {
    row_split_pass <- ann_data[, split_heatmap]
  }
  if (!is.null(annotationData)) {
    if (length(colList) == 0) {
      colorSetNum <- 1
      for (annot in annotationColNames) {
        if (is.numeric(annotationData[, annot])) {
          t1min <- min(annotationData[, annot], na.rm = TRUE)
          t1max <- max(annotationData[, annot], na.rm = TRUE)
          colList[[annot]] <- circlize::colorRamp2(c(t1min, t1max),
                                                   c("white", "blue"))
        } else {
          if (is.factor(annotationData[, annot][!is.na(annotationData[, annot])])) {
            condLevels <- levels(annotationData[, annot][!is.na(annotationData[, annot])])
          } else {
            condLevels <- unique(annotationData[, annot][!is.na(annotationData[, annot])])
          }
          if (length(condLevels) > 8) {
            colors <- distinctColors(length(condLevels))
          } else {
            colors <- RColorBrewer::brewer.pal(8, colorSets[colorSetNum])
            colorSetNum <- colorSetNum + 1
          }
          colList[[annot]] <- stats::setNames(colors[seq_along(condLevels)],
                                              condLevels)
        }
      }
    }
    topha2 <- ComplexHeatmap::HeatmapAnnotation(
      df = data.frame(annotationData),
      col = colList, height = grid::unit(0.4 * length(annotationColNames), "cm"),
      show_legend = TRUE, show_annotation_name = TRUE)
    return(ComplexHeatmap::draw(
      ComplexHeatmap::Heatmap(sigresults, column_title = plot_title,
                              name = name,
                              show_column_names = showColumnNames,
                              col = choose_color,
                              show_row_names = showRowNames,
                              top_annotation = topha2,
                              row_split = row_split_pass,
                              column_order = column_order,
                              cluster_columns = TRUE,
                              ...),
      annotation_legend_side = "bottom"))
  }
}