#' Parse a MultiExperimentAssay object and extract the elements as data.frames
#'
#' This function takes an animalcules-formatted \code{MultiAssayExperiment}
#' object and parses it to extract a named assay alongside the taxonomy and
#' metadata.
#'
#' @param dat Either a \code{MultiAssayExperiment} or a
#'   \code{SummarizedExperiment} object.
#' @param which_assay Character string. This refers to the assay to be extracted
#'   from the MultiAssayExperiment object if \code{type = "MAE"}. Does not need
#'   to be specified if \code{type = "SE"}. Default is \code{NULL}.
#' @param type One of "MAE" denoting a \code{MultiAssayExperiment} or "SE"
#'   denoting a \code{SummarizedExperiment}.
#'
#' @returns Returns a list of 3 named data.frame elements, `counts`, `sam`, and
#'   `tax` denoting the counts data, sample metadata table, and taxonomy table,
#'   respectively.
#'
#' @export
#' @importFrom rlang .data
#' @import MultiAssayExperiment
#'
#' @examples
#' in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
#' out <- parse_MAE_SE(in_dat)
#' head(out$tax)
#' head(out$sam)
#' head(out$counts)
#' 
#' out2 <- parse_MAE_SE(in_dat[["MicrobeGenetics"]],
#'                      which_assay = "MGX", type = "SE")
#' 

parse_MAE_SE <- function(dat, which_assay = NULL, type = "MAE") {
  if (type == "MAE") {
    if (is.null(which_assay)) {
      which_assay <- names(MultiAssayExperiment::assays(dat))[1]
    }
    if (!methods::is(dat, "MultiAssayExperiment")) {
      stop("Input must be a MultiAssayExperiment")
    }
    microbe <- dat[[which_assay]]
  } else if (type == "SE") {
    microbe <- dat
  } else stop("type must be one of 'SE' or 'MAE'")
  
  # Extract metadata, taxonomic info, and counts
  tax_table <- as.data.frame(SummarizedExperiment::rowData(microbe))
  sam_table <- as.data.frame(SummarizedExperiment::colData(microbe))
  counts_table <- as.data.frame(SummarizedExperiment::assay(
    microbe, "MGX"))[, rownames(sam_table)]
  
  list(counts = counts_table,
       sam = sam_table,
       tax = tax_table) %>%
    return()
}
