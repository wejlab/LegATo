#' Create a table summarizing reads aggregated by grouping variables on a unit
#' 
#' This function takes a \code{MultiAssayExperiment} of microbial read counts
#' and aggregates them by one or more grouping vars within a unit.
#' 
#' @inheritParams plot_stacked_bar
#' @param group_vars A character string or character vector of covariates
#' found in \code{colData(dat)} to use in grouping counts. The variables
#' should be listed in order of desired grouping.
#' 
#' @return A \code{data.frame} of the grouping columns, mean_reads, sd_reads,
#' min_reads, max_reads and num_total.
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
#' out <- get_summary_table(in_dat, c("Group", "Subject"))
#' head(out)
#'

get_summary_table <- function(dat, group_vars) {
  microbe <- parse_MAE_SE(dat)
  counts_table <- microbe$counts
  tax_table <- microbe$tax
  sam_table <- microbe$sam
  
  reads_tbl <- counts_table %>%
    tidyr::as_tibble() %>%
    dplyr::mutate(species = tax_table$species) %>%
    dplyr::relocate(.data$`species`) %>%
    tidyr::pivot_longer(cols = tidyr::all_of(colnames(counts_table)),
                        names_to = "Sample_id",
                        values_to = "Abundance") %>%
    dplyr::group_by(.data$Sample_id) %>%
    dplyr::summarise(TotalAbundance = sum(.data$Abundance)) %>%
    dplyr::left_join(tibble::rownames_to_column(sam_table, "Sample_id"),
                     by = "Sample_id") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_vars))))
  
  summary_tab <- reads_tbl %>%
    dplyr::summarise(mean_reads = base::mean(.data$TotalAbundance),
                     sd_reads = stats::sd(.data$TotalAbundance),
                     min_reads = base::min(.data$TotalAbundance),
                     max_reads = base::max(.data$TotalAbundance),
                     num_total = dplyr::n()) %>%
    as.data.frame()
  
  return(summary_tab)
}
