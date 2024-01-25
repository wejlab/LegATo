#' Clean up taxon names in a MultiAssayExperiment
#'
#' This functional is an optional method for fixing up taxon names in a
#' \code{MultiAssayExperiment} to be run before \code{filter_animalcules_MAE}.
#' Specifically it removes brackets from species names, replaces species labeled
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
#' clean_animalcules_MAE(in_dat)
#' 

clean_animalcules_MAE <- function(dat) {
  parsed <- parse_MAE_SE(dat, which_assay = "MicrobeGenetics", type = "MAE")
  tax_table <- parsed$tax
  sam_table <- parsed$sam
  counts_table <- parsed$counts

  se_rowData  <- tax_table %>%
    dplyr::mutate("species" = stringr::str_remove_all(.data$species, "\\[|\\]"))
  
  # Replace species that are others with "sp."
  ind <- se_rowData$species == "others"
  se_rowData$species[ind] <- paste(se_rowData$genus[ind], "sp.", sep = " ")
  # Replace genus that are others with species information
  ind <- se_rowData$genus == "others" & !is.na(se_rowData$genus)
  all_split <- sapply(strsplit(se_rowData$species, " "), function(x) x[[1]])
  se_rowData$genus[ind] <- all_split[ind]
  
  #ISSUE: DUplicates because we need to consolidate down to species level
  
  # Fix rownames
  if (all.equal(rownames(tax_table), rownames(counts_table))) {
    to_rownames <- stringr::str_replace_all(se_rowData$species, "_", " ")
    rownames(se_rowData) <- rownames(counts_table) <- to_rownames
    # Alphabetize
    se_rowData_final <- se_rowData[order(to_rownames),]
    counts_table_final <- counts_table[order(to_rownames),]
  } else (stop("Mismatch of rownames in datasets"))

  se_mgx <- counts_table %>% base::data.matrix() %>% S4Vectors::SimpleList() %>%
    magrittr::set_names("MGX")

  microbe_se <- SummarizedExperiment::SummarizedExperiment(
    assays = se_mgx, colData = sam_table, rowData = se_rowData)
  MAE_out <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList(MicrobeGenetics = microbe_se),
    colData = sam_table)
  
  return(MAE_out)
}

