#' Obtain a data.frame of ordered taxa abundances at a given level
#'
#' This function takes a \code{MultiAssayExperiment} object and returns a
#' data.frame of the present taxa at a user-supplied taxonomy level,
#' and outputs the average abundances of the taxa. 
#'
#' @inheritParams plot_stacked_bar
#'
#' @return A \code{data.frame}
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
#' out <- get_top_taxa(in_dat, "genus")
#' out
#'

get_top_taxa <- function(dat, taxon_level = "genus") {
  microbe <- parse_MAE_SE(dat)
  all_relabu <- microbe$counts |>
    animalcules::upsample_counts(microbe$tax, taxon_level) %>%
    animalcules::counts_to_relabu() %>%
    tibble::rownames_to_column(var = "taxon") |>
    dplyr::rowwise("taxon") |>
    # Sum everything but the first columm (taxon name)
    dplyr::summarise(
      allmeans = mean(
        dplyr::c_across(dplyr::all_of(colnames(microbe$counts)))),
      .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$allmeans)) |>
    as.data.frame()
  
  return(all_relabu)
}
