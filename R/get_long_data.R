#' Create a long data.frame from a MultiAssayExperiment counts object
#'
#' This function takes a \code{MultiAssayExperiment} object and a specified
#' taxon level of interest and creates a long \code{data.frame} that can be used
#' more easily for plotting counts data.
#'
#' @inheritParams plot_stacked_bar
#' @param log logical. Indicate whether an assay returned should be the log of
#'   whichever assay is specified in \code{"output_name"}. If
#'   \code{counts_to_CPM = TRUE} as well, then a log CPM assay will also be
#'   created. Default is \code{FALSE}.
#' @param counts_to_CPM logical. This argument only applies if the
#'   \code{input_type} is a counts assay. If \code{TRUE}, then the output assays
#'   will include a normalized CPM assay. If \code{log = TRUE} as well, then a
#'   log CPM assay will also be created. Default is \code{TRUE}.
#'
#' @return A \code{data.frame} consisting of the counts data, taxa, and
#'   metadata.
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
#' out <- get_long_data(in_dat, "genus", log = TRUE, counts_to_CPM = TRUE)
#' head(out)
#' 

get_long_data <- function(dat, taxon_level, log = FALSE,
                        counts_to_CPM = FALSE) {
  SE_obj <- dat[["MicrobeGenetics"]]
  assay_name <- names(SummarizedExperiment::assays(SE_obj))[1]
  if (log | counts_to_CPM) {
    SE_obj <- TBSignatureProfiler::mkAssay(SE_obj,
                                           input_name = assay_name,
                                           output_name = "assay",
                                           log = log, counts_to_CPM = counts_to_CPM)
    if (log) which_assay = "log_assay"
    if (counts_to_CPM) which_assay = "assay_CPM"
    if (log && counts_to_CPM) which_assay = "log_assay_CPM"
  }
  
  microbe <- parse_MAE_SE(SE_obj, which_assay = which_assay, type = "SE")
  
  output <- animalcules::upsample_counts(microbe$counts, microbe$tax, taxon_level) %>%
    tibble::rownames_to_column("taxon") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(colnames(microbe$counts)),
                        names_to = "Sample_id", values_to = "Abundance") %>%
    dplyr::left_join(tibble::rownames_to_column(microbe$sam, "Sample_id"),
                     by = "Sample_id") %>%
    dplyr::relocate("taxon", colnames(microbe$sam)) %>%
    as.data.frame()
  return(output)
}
