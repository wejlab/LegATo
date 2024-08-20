#' Clean up taxon names in a MultiAssayExperiment
#'
#' This functional is an optional method for fixing up taxon names in a
#' \code{MultiAssayExperiment} to be run before \code{filter_MAE}.
#' Specifically, it removes brackets from species names, replaces species labeled
#' as "others" with "sp." and finally replaces underscores with spaces.
#'
#' @inheritParams plot_stacked_bar
#'
#' @returns An animalcules-formatted \code{MultiAssayExperiment} object with
#'   cleaned-up taxonomy nomenclature.
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
#' clean_MAE(in_dat)
#'

clean_MAE <- function(dat) {
  # Extract data
  parsed <- parse_MAE_SE(dat, which_assay = "MicrobeGenetics", type = "MAE")
  tax_table <- parsed$tax
  counts_table <- parsed$counts
  # Preliminary fixing of species names
  if(!all(c("genus", "species") %in% colnames(tax_table))) {
    stop("Columns 'genus' and 'species' must be present in taxonomy table.")
  } else if(!all.equal(rownames(tax_table), rownames(counts_table))) {
    stop("MAE could not be parsed correctly. Please ensure correct formatting.")
  }
  # Only need to alter species names in tax_table initially
  se_rowData  <- tax_table %>%
    dplyr::mutate("species" = stringr::str_remove_all(.data$species, "\\[|\\]"),
                  "species" = stringr::str_replace_all(.data$species, "_", " "))
  ## Replace species that are others with "sp."
  ind <- se_rowData$species == "others"
  se_rowData$species[ind] <- paste(se_rowData$genus[ind], "sp.", sep = " ")
  ## Replace genus that are others with species information
  ind <- se_rowData$genus == "others" & !is.na(se_rowData$genus)
  all_split <- plyr::aaply(strsplit(se_rowData$species, " "), 1,
                           function(x) x[[1]])
  se_rowData$genus[ind] <- all_split[ind]
  # Consolidate tax and counts to species level
  up_counts <- animalcules::upsample_counts(counts_table, se_rowData,
                                            higher_level = "species")
  up_rowData <- se_rowData %>%
    tibble::remove_rownames() %>%
    dplyr::distinct(.data$species, .keep_all = TRUE)
  ind <- match(rownames(up_counts), up_rowData$species)
  up_rowData <- up_rowData[ind, ]
  rownames(up_rowData) <- up_rowData$species
  
  MAE_out <- create_formatted_MAE(counts_dat = up_counts,
                                  tax_dat = up_rowData,
                                  metadata_dat = parsed$sam)
  return(MAE_out)
}

