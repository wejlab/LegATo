plot_spaghetti <- function(dat,
                           taxon_level,
                           covariate_time,
                           covariate_1,
                           unit_var,
                           which_taxon,
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
                  subtitle = subtitle,
                  ylabel = "Relative Abundance (log CPM)") +
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
