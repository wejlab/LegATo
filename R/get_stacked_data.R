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
