#' Parse a MultiExperimentAssay object and extract the elements as data.frames
#'
#' This function takes an animalcules-formatted \code{MultiAssayExperiment}
#' object and parses it to extract a named assay alongside the taxonomy and
#' metadata.
#'
#' @param dat Either a \code{MultiAssayExperiment} or a
#'   \code{SummarizedExperiment} object.
#' @param which_experiment Character string. If \code{type = "MAE"}, then
#' this is the name of the experiment to be accessed. If \code{NULL},
#' will default to the first available experiment.
#' @param which_assay Character string. Regardless of whether 
#' \code{type = "MAE"} or \code{"SE"}, this is the name of the selected
#' \code{SummarizedExperiment} object. If \code{NULL}, defaults to first listed.
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

parse_MAE_SE <- function(dat, which_experiment = NULL,
                         which_assay = NULL, type = "MAE") {
    if (type == "MAE") {
        if (is.null(which_experiment)) {
            which_experiment <- names(MultiAssayExperiment::assays(dat))[1]
        } 
        if (!methods::is(dat, "MultiAssayExperiment")) {
            stop("Input must be a MultiAssayExperiment")
        }
        microbe <- dat[[which_experiment]]
    } else if (type == "SE") {
        microbe <- dat
    } else stop("type must be one of 'SE' or 'MAE'")
    
    if (is.null(which_assay)) {
        which_assay <- names(SummarizedExperiment::assays(microbe))[1]
    }
    
    # Extract metadata, taxonomic info, and counts
    tax_table <- as.data.frame(SummarizedExperiment::rowData(microbe))
    sam_table <- as.data.frame(SummarizedExperiment::colData(microbe))
    
    SE_assay <- SummarizedExperiment::assay(microbe, which_assay)
    counts_table <- as.data.frame(SE_assay)[, rownames(sam_table)]
    
    list(counts = counts_table,
         sam = sam_table,
         tax = tax_table) %>%
        return()
}
