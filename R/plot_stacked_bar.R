#' Plot a stacked bar chart of microbial relative abundance
#'
#' This function takes a \code{MultiAssayExperiment} object and returns a
#' stacked bar plot of microbe relative abundances. The function takes a single
#' covariate as an optional variable by which to create multiple gridded plots.
#'
#' If further manipulation of specific parameters is desired, users can add
#' \code{ggplot2} function calls to the output of the function.
#'
#' @param dat A \code{MultiAssayExperiment} object specially formatted as an
#'   animalcules output.
#' @param taxon_level Character string indicating the level of taxonomy to
#'   aggregate the counts data. Must be the name of a column in
#'   \code{MultiAssayExperiment::rowData(dat)}.
#' @param covariate_1 Character string giving the name of a column in
#'   \code{MultiAssayExperiment::colData(dat)} on which to create multiple
#'   plots. The default is \code{NULL}.
#' @param covariate_time Character string giving the name of the discrete
#'   time-based covariate in the metadata to group the points on the bar chart.
#' @param palette_input A character vector of colors that is at minimum the same
#'   length of the number of taxa (specified with \code{taxon_level}).
#'   The default is \code{NULL} and relies on \code{ggplot2}'s default scheme.
#' @param title A character string providing the plot title.
#' @param subtitle A character string providing the plot subtitle. The default
#'   is \code{NULL}.
#'
#' @return A \code{ggplot2} plot.
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' print("My example here")
#'

plot_stacked_bar <- function(dat,
                             taxon_level,
                             covariate_1 = NULL,
                             covariate_time,
                             palette_input= NULL,
                             title = paste("Relative abundance at", taxon_level, "level"),
                             subtitle = NULL) {
  relabu_table <- get_stacked_data(dat, taxon_level, covariate_1, covariate_time)
  
  taxa_ordered <- get_top_taxa(dat, taxon_level) %>%
    dplyr::pull(.data$taxa)
  
  myplot <- relabu_table %>%
    dplyr::mutate("Taxon" = factor(.data$taxon, levels = taxa_ordered),
                  "Timepoint" = factor(covariate_t)) %>%
    dplyr::rename("Relative abundance" = .data$value) %>%
    ggplot2::ggplot(ggplot2::aes(fill = .data$Taxon, x = .data$Timepoint, y = .data$`Relative abundance`)) + 
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = title,
                  subtitle = subtitle) +
    ggplot2::theme(legend.position = "bottom",
                   #axis.title.x = element_blank(), axis.text.x = element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank()) +
    ggplot2::coord_flip()
  
  if (!is.null(palette_input)) {
    myplot <- myplot +
      ggplot2::scale_fill_manual(values = palette_input)
  }
  if (!is.null(covariate_1)) {
    myplot <- myplot +
    ggplot2::facet_grid(~.data$covariate_1)
  }
  
  return(myplot)
}
