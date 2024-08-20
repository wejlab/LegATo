
test_models_lmm <- function(tn, input_df, unit_var, fixed_cov,
                            plotsave_loc, plot_terms, plot_out, ...) {
  filt_df <- input_df %>% 
    dplyr::mutate("unit_var" = as.factor(input_df[, unit_var])) %>%
    dplyr::filter(.data$taxon == tn) %>%
    dplyr::arrange(.data$unit_var)
  group_vars <- paste("Abundance ~", paste(fixed_cov, collapse = " + ")) %>%
    paste0(" + (1|", unit_var, ")") %>%
    stats::formula()
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    message("The 'lmerTest' package is not installed",
            "Please install it from CRAN to use this function.")
  }
  complex <- lmerTest::lmer(group_vars, data = filt_df,
                        na.action = stats::na.omit)
  if (requireNamespace("broom", quietly = TRUE) &&
      requireNamespace("broom.mixed", quietly = TRUE)) {
    sum_comp <- broom.mixed::tidy(complex) %>%
      dplyr::filter(is.na(.data$group)) %>%
      dplyr::select(-"group")
  } else {
    message("The 'broom' and/or 'broom.mixed' packages are not installed",
            "Please install them to use this function.")
  }
  
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


#' Compute linear mixed-effects models (LMM) on longitudinal microbiome data
#'
#' This function takes an animalcules-formatted \code{MultiAssayExperiment} and
#' runs an independent LMM model for each taxon. The model predicts taxon log
#' CPM abundance as a product of fixed-effects covariates with a random effect,
#' usually the unit on which repeated measurements were taken. Note, the 'broom',
#' 'lmerTest', and 'broom.mixed' packages are required to use this function; they can
#' be downloaded from CRAN.
#'
#' P-values are adjusted for the model coefficients within each taxon. The
#' following methods are permitted: \code{c("holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none")}
#' 
#' 
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
#'   filter_MAE()
#' out <- run_lmm_model(dat, taxon_level = "genus", unit_var = "Subject",
#'                      fixed_cov = c("HIVStatus", "timepoint"))
#' head(out)
#'

run_lmm_model <- function(dat,
                          taxon_level = "genus",
                          unit_var,
                          fixed_cov,
                          p_adj_method = "fdr",
                          plot_out = FALSE,
                          plotsave_loc = ".",
                          plot_terms = NULL,
                          ...) {
  input_df <- get_long_data(dat, taxon_level, log = TRUE, counts_to_CPM = TRUE)
  all_tn <- get_top_taxa(dat, taxon_level) %>% dplyr::pull("taxon")
  n <- length(all_tn)
  storage <- plyr::llply(all_tn, test_models_lmm, input_df = input_df,
                         fixed_cov = fixed_cov,
                         unit_var = unit_var,
                         plot_out = plot_out,
                         plotsave_loc = plotsave_loc,
                         plot_terms = plot_terms,
                         ...) %>%
    data.table::rbindlist() %>%
    dplyr::arrange(.data$Coefficient) %>%
    dplyr::group_by(.data$Coefficient) %>%
    dplyr::mutate("Adj p-value" = stats::p.adjust(
      .data$`Unadj p-value`, method = p_adj_method)) %>%
    as.data.frame()
  return(storage)
}
