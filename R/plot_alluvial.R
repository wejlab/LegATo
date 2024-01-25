#' Plot an alluvial diagram of microbial relative abundance
#'
#' This function takes a \code{MultiAssayExperiment} object and returns a
#' alluvial diagram of microbe relative abundances. The function takes a single
#' covariate as an optional variable by which to create a grid of multiple
#' plots.
#'
#' If further manipulation of specific parameters is desired, users can add
#' \code{ggplot2} function calls to the output of the function.
#'
#' @inheritParams plot_stacked_bar
#'
#' @return A \code{ggplot2} plot.
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") %>% readRDS()
#' plot_alluvial(in_dat, taxon_level = "family", covariate_1 = "Group", covariate_time = "Month",
#'               palette_input = rainbow(25))
#'

plot_alluvial <- function(dat,
                          taxon_level,
                          covariate_1 = NULL,
                          covariate_time,
                          palette_input = NULL,
                          title = paste("Relative abundance at", taxon_level, "level"),
                          subtitle = NULL) {
  relabu_table <- get_stacked_data(dat, taxon_level, covariate_1, covariate_time)
  
  taxa_ordered <- get_top_taxa(dat, taxon_level) %>%
    dplyr::pull(.data$taxon)
  
  myplot <- relabu_table %>%
    dplyr::mutate("Taxon" = factor(.data$taxon, levels = taxa_ordered),
                  "Timepoint" = as.factor(.data$covariate_t),
                  "covariate_1" = factor(.data$covariate_1)) %>%
    dplyr::rename("Relative abundance" = .data$value) %>%
    dplyr::select("covariate_1", "Timepoint", "Taxon", "Relative abundance") %>%
    ggplot2::ggplot(ggplot2::aes(y = .data$`Relative abundance`,
                                 x = .data$Timepoint, alluvium = .data$Taxon)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = .data$Taxon,
                                           color = .data$Taxon),
                              width = 1/4, alpha = 0.7, decreasing = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Timepoint", y = "Relative abundance",
                  title = title, subtitle = subtitle) +
    ggplot2::theme(legend.position = "bottom",
                   legend.title = ggplot2::element_blank())
  if (!is.null(palette_input)) {
    myplot <- myplot +
      ggplot2::scale_fill_manual(values = palette_input) +
      ggplot2::scale_color_manual(values = palette_input)
  }
  if (!is.null(covariate_1)) {
    myplot <- myplot +
      ggplot2::facet_grid(~.data$covariate_1)
  }
  
  return(myplot)
}

