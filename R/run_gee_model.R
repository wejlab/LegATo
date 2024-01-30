# If emmeans doesn't work, may need to step through function
mk_gee_plot <- function(this_coef, complex, tn, unit_var,
                        plotsave_loc, plot_terms, ...) {
  term_input <- unique(c(this_coef, plot_terms))
  p <- plot(ggeffects::ggemmeans(complex, terms = term_input)) + 
    ggplot2::labs(subtitle = "Estimated marginal means",
                  title = stringr::str_replace(tn, "_", " ")) +
    ggplot2::ylab("Taxon abundance (log CPM)") +
    ggplot2::theme_classic()
  filename <- file.path(plotsave_loc, paste0(tn, "_", this_coef, ".png"))
  ggplot2::ggsave(filename = filename, p, ...)
}

test_models_gee <- function(tn, input_df, unit_var, fixed_cov,
                            corstr, plot_out,
                            plotsave_loc, plot_terms, ...) {
  filt_df <- input_df %>% 
    dplyr::mutate("unit_var" = as.factor(input_df[, unit_var])) %>%
    dplyr::filter(.data$taxon == tn) %>%
    dplyr::arrange(.data$unit_var)
  group_vars <- paste("Abundance ~", paste(fixed_cov, collapse = " + ")) %>%
    stats::formula()
  complex <- geepack::geeglm(group_vars, id = unit_var, data = filt_df,
                             na.action = stats::na.omit, family = "gaussian",
                             corstr = corstr)
  sum_comp <- broom::tidy(complex, conf.int = TRUE)

  if (plot_out) {
    plyr::a_ply(fixed_cov, 1, mk_gee_plot, complex = complex, tn = tn,
                plotsave_loc = plotsave_loc, plot_terms = plot_terms, ...)
  }
  res_out <- sum_comp %>%
    dplyr::mutate("Taxon" = tn) %>%
    dplyr::rename("Coefficient" = "term",
                  "Coefficient Estimate" = "estimate",
                  "Lower 95% CI" = "conf.low",
                  "Upper 95% CI" = "conf.high",
                  "Standard Error" = "std.error",
                  "Statistic" = "statistic",
                  "Pr(>|W|)" = "p.value") %>%
    as.data.frame()
  return(res_out)
}

#' Compute Generalized Estimating Equations (GEEs) on longitudinal microbiome
#' data
#'
#' This function takes an animalcules-formatted \code{MultiAssayExperiment} and
#' runs an independent GEE model for each taxon. The model predicts taxon log
#' CPM abundance as a product of fixed-effects covariates conditional on a
#' grouping ID variable, usually the unit on which repeated measurements were
#' taken. This modeling approach works best with small datasets that multiple
#' samples across many (>40) clusters/units.
#'
#' P-values are adjusted for the model coefficients within each taxon. The
#' following methods are permitted: \code{c("holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none")}
#'
#' @inheritParams test_hotelling_t2
#' @param fixed_cov A character vector naming covariates to be tested.
#' @param corstr A character string specifying the correlation structure. The
#'   following are permitted: '"independence"', '"exchangeable"', '"ar1"',
#'   '"unstructured"'.
#' @param p_adj_method A character string specifying the correction method. Can
#'   be abbreviated. See details. Default is \code{"fdr"}.
#' @param plot_out Logical indicating whether plots should be output alongside
#'   the model results. Default is \code{FALSE}.
#' @param plotsave_loc A character string giving the folder path to save plot
#'   outputs. This defaults to the current working directory.
#' @param plot_terms Character vector. Which terms should be examined in the
#'   plot output? Can overlap with the \code{fixed_cov} inputs.
#' @param ... Further arguments passed to \code{ggsave} for plot creation.
#'
#' @export
#' @importFrom rlang .data
#'
#' @returns A \code{data.frame} of modeling results.
#' @examples
#' in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |>
#'               readRDS()
#' out <- run_gee_model(in_dat, taxon_level = "genus", unit_var = "Subject",
#'                      fixed_cov = c("HairLength", "Age", "Group", "Sex"),
#'                      corstr = "ar1")
#' head(out)
#' 

run_gee_model <- function(dat,
                          taxon_level = "genus",
                          unit_var,
                          fixed_cov,
                          corstr = "ar1",
                          p_adj_method = "fdr",
                          plot_out = FALSE,
                          plotsave_loc = ".",
                          plot_terms = NULL,
                          ...) {
  input_df <- get_long_data(dat, taxon_level, log = TRUE, counts_to_CPM = TRUE)
  all_tn <- get_top_taxa(dat, taxon_level) %>% dplyr::pull("taxon")
  n <- length(all_tn)
  storage <- plyr::llply(all_tn, test_models_gee, input_df = input_df,
                         fixed_cov = fixed_cov,
                         unit_var = unit_var,
                         corstr = corstr,
                         plot_out = plot_out,
                         plotsave_loc = plotsave_loc,
                         plot_terms = plot_terms, ...) %>%
    data.table::rbindlist() %>%
    dplyr::arrange(.data$Coefficient) %>%
    dplyr::group_by(.data$Coefficient) %>%
    dplyr::mutate("Adj p-value" = stats::p.adjust(.data$`Pr(>|W|)`, method = p_adj_method)) %>%
    dplyr::rename("Unadj p-value" = .data$`Pr(>|W|)`) %>%
    as.data.frame()
  return(storage)
}