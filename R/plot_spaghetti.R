#' Plot a spaghetti volatility plot of microbial abundance for a given taxon
#'
#' This function takes a \code{MultiAssayExperiment} object and returns a
#' spaghetti plot of microbial abundance delineated by a unit, such as
#' a subject.
#'
#' If further manipulation of specific parameters is desired, users can add
#' \code{ggplot2} function calls to the output of the function.
#'
#' @inheritParams plot_stacked_bar
#' @param covariate_1 Character string, the name of the covariate in `dat`
#' by which to color and group samples. Default is \code{NULL}.
#' @param unit_var Character string, the name of the column delineating the
#' unit on which the microbial abundances are changing over time. This is
#' likely something akin to a subject that repeated measurements are made on.
#' @param which_taxon Character string, the name of the taxon to plot at the
#' specified taxon level.
#'
#' @return A \code{ggplot2} plot.
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
#' all_taxa <- get_top_taxa(in_dat, "phylum")
#' plot_spaghetti(in_dat, taxon_level = "phylum", covariate_1 = "Group", covariate_time = "Month",
#'               unit_var = "Subject", which_taxon = all_taxa$taxon[1],
#'               palette_input = rainbow(25))
#'

plot_spaghetti <- function(dat,
                           covariate_time,
                           covariate_1 = NULL,
                           unit_var,
                           taxon_level,
                           which_taxon,
                           palette_input= NULL,
                           title = "Spaghetti Plot",
                           subtitle = NULL) {
  input_data <- get_long_data(dat, taxon_level, log = TRUE,
                              counts_to_CPM = TRUE) %>%
    dplyr::filter(.data$taxon == which_taxon)
  p <- input_data %>%
    ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym(covariate_time),
                                 y = .data$Abundance,
                                 group = !!rlang::sym(unit_var),
                                 color = !!rlang::sym(covariate_1))) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_point(alpha = 0.7, size = 0.5) +
    ggplot2::labs(title = "Spaghetti plot of microbe abundance",
                  subtitle = which_taxon) +
    ggplot2::theme_bw()
  
  if (!is.null(palette_input)) {
    p <- p + ggplot2::scale_color_manual(values = palette_input)
  }
  if(!is.null(covariate_1)) {
    p <- p + ggplot2::facet_wrap(rlang::sym(covariate_1))
  }
  return(p)
}
