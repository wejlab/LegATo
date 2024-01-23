#' Plot a stacked area chart of microbial relative abundance
#'
#' This function takes a \code{MultiAssayExperiment} object and returns a
#' stacked area chart of microbe relative abundances. The function takes a single
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
#' print("My example here")
#' 

plot_stacked_area <- function(dat,
                              taxon_level,
                              covariate_1 = NULL,
                              covariate_time,
                              palette_input = NULL,
                              title = paste("Relative abundance at", taxon_level, "level"),
                              subtitle = NULL) {
  
  relabu_table <- get_stacked_data(dat, taxon_level,
                                   covariate_1, covariate_time)
  
  taxa_ordered <- get_top_taxa(dat, taxon_level) %>%
    dplyr::pull(taxa)
  
  p <- relabu_table %>%
    dplyr::mutate("Taxon" = factor(.data$taxon, levels = taxa_ordered),
                  "Time point" = .data$covariate_t) %>%
    dplyr::rename("Relative abundance" = .data$value) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$`Time point`,
                                 y = .data$`Relative abundance`,
                                 fill = .data$`Taxon`)) + 
    ggplot2::geom_area(alpha = 0.7, linewidth = .5, colour = "white") +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Time point", y = "Relative abundance",
                  title = title,
                  subtitle = subtitle) +
    ggplot2::theme(legend.position = "bottom",
                   legend.title = ggplot2::element_blank())
  
  if(!is.null(covariate_1)) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(covariate_1))
  } 
  if (!is.null(palette_input)) {
    p <- p + ggplot2::scale_fill_manual(values = palette_input)
  }
  return(p)
}
