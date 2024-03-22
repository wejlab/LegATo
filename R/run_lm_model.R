test_models_lm <- function(tn, input_df, fixed_cov,
                            plotsave_loc, plot_terms, plot_out, ...) {
  filt_df <- input_df %>%
    dplyr::filter(.data$taxon == tn)
  group_vars <- paste("Abundance ~", paste(fixed_cov, collapse = " + ")) %>%
    stats::formula()
  complex <- stats::lm(group_vars, data = filt_df,
                       na.action = stats::na.omit)
  #confidence <- lme4::confint.merMod(complex, quiet = FALSE) %>%
  #  as.data.frame() %>%
  #  tibble::rownames_to_column("term")
  sum_comp <- broom.mixed::tidy(complex)
    # dplyr::left_join(confidence, by = "term")
  
  if (plot_out) {
    plyr::a_ply(fixed_cov, 1, mk_gee_plot, complex = complex, tn = tn,
                plotsave_loc = plotsave_loc, plot_terms = plot_terms, ...)
  }
  res_out <- sum_comp %>%
    dplyr::mutate("Taxon" = tn) %>%
    dplyr::rename("Coefficient" = "term",
                  "Coefficient Estimate" = "estimate",
                  #              "Lower 95% CI" = "2.5 %",
                  #              "Upper 95% CI" = "97.5 %",
                  "Standard Error" = "std.error",
                  "Statistic" = "statistic",
                  "Unadj p-value" = "p.value") %>%
    as.data.frame()
  return(res_out)
}

#' Compute linear models (LMM) on microbiome data
#'
#' This function takes an animalcules-formatted \code{MultiAssayExperiment} and
#' runs an independent linear model for each taxon. The model predicts taxon log
#' CPM abundance as a product of user-specified covariates. This model can
#' be used for general microbiome analyses without repeated measures data.
#'
#' P-values are adjusted for the model coefficients within each taxon. The
#' following methods are permitted: \code{c("holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none")}
#'
#' @inheritParams run_gee_model
#'
#' @export
#' @importFrom rlang .data
#'
#' @returns A \code{data.frame} of modeling results.
#' @examples
#' dat <- system.file("extdata/MAE.RDS", package = "LegATo") |>
#'   readRDS() |>
#'   filter_animalcules_MAE(0.05)
#' out <- run_lm_model(dat, fixed_cov = c("timepoint", "HIVStatus"),
#'                     plot_out = FALSE)
#' head(out)
#'
run_lm_model <- function(dat,
                          taxon_level = "genus",
                          fixed_cov,
                          p_adj_method = "fdr",
                          plot_out = FALSE,
                          plotsave_loc = ".",
                          plot_terms = NULL,
                          ...) {
  input_df <- get_long_data(dat, taxon_level, log = TRUE, counts_to_CPM = TRUE)
  all_tn <- get_top_taxa(dat, taxon_level) %>% dplyr::pull("taxon")
  n <- length(all_tn)
  storage <- plyr::llply(all_tn, test_models_lm, 
                         input_df = input_df,
                         fixed_cov = fixed_cov,
                         plot_out = plot_out,
                         plotsave_loc = plotsave_loc,
                         plot_terms = plot_terms) %>%
    data.table::rbindlist() %>%
    dplyr::arrange(.data$Coefficient) %>%
    dplyr::group_by(.data$Coefficient) %>%
    dplyr::mutate("Adj p-value" = stats::p.adjust(
      .data$`Unadj p-value`, method = p_adj_method)) %>%
    as.data.frame()
  return(storage)
}
