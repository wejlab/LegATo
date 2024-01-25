#' Create a long data.frame with grouped abundances from a MultiAssayExperiment counts object
#' 
#' This function takes a \code{MultiAssayExperiment} object and a specified
#' taxon level of interest and creates a long \code{data.frame} that can be used
#' more easily for plotting counts data in a stacked bar plot or a stacked area
#' chart. The function groups taxa and computes relative abundance within taxa strata.
#' 
#' @inheritParams plot_spaghetti
#' 
#' @return A \code{data.frame} consisting of the counts data, taxa, and metadata.
#'
#' @export
#' @importFrom rlang .data
#' 
#' @examples
#' in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
#' get_stacked_data(in_dat, "genus", covariate_1 = "Sex", covariate_time = "Month")
#'

get_stacked_data <- function(dat, taxon_level = "genus",
                             covariate_1,
                             covariate_time) {
  microbe <- parse_MAE_SE(dat)
  
  if (is.null(covariate_1)) {
    relabu_table <- microbe$counts %>%
      animalcules::upsample_counts(., microbe$tax, taxon_level) %>%
      animalcules::counts_to_relabu() %>% t() %>% as.data.frame() %>%
      # Add grouping vars
      dplyr::mutate("covariate_t" = microbe$sam[, covariate_time]) %>%
      tidyr::pivot_longer(!"covariate_t", names_to = "taxon") %>%
      S4Vectors::aggregate(. ~ taxon + `covariate_t`, ., mean)
  } else {
    relabu_table <- microbe$counts %>%
      animalcules::upsample_counts(., microbe$tax, taxon_level) %>%
      animalcules::counts_to_relabu() %>% t() %>% as.data.frame() %>%
      # Add grouping vars
      dplyr::mutate("covariate_1" = microbe$sam[, covariate_1],
                    "covariate_t" = microbe$sam[, covariate_time]) %>%
      tidyr::pivot_longer(!c("covariate_1", "covariate_t"), names_to = "taxon") %>%
      S4Vectors::aggregate(. ~ taxon + covariate_1 + covariate_t, ., mean)
  }
  # formula for covariates
  return(relabu_table)
}
