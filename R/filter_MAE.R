
#' Filter a MultiAssayExperiment object to keep a top percentage of taxa
#' 
#' This function takes an animalcules-formatted \code{MultiAssayExperiment} (MAE)
#' object and identifies all taxa at the OTU level of choice that exhibit a
#' relative abundance greater than or equal to a relative abundance percent
#' threshold, \code{relabu_threshold}, in at least \code{occur_pct_cutoff}\% of
#' the total samples. After filtration, taxa across the specified OTU level and
#' all downstream levels are then consolidated into the category "Other".
#' @inheritParams plot_stacked_bar
#' @param relabu_threshold A double(percentage) between 0 and 100, representing
#' the relative abundance criterion that all OTUs should meet to be retained.
#' The smaller the threshold, the fewer the OTUs will be retained. Default is
#' 3\%.
#' @param occur_pct_cutoff A double (percentage) between 0 and 100
#' representing the percent cutoff for how many OTUs must meet the
#' \code{relabu_threshold} across the samples to be retained. It is wise to
#' keep the number of samples in mind when setting this parameter. Default is 5\%.
#' @returns An animalcules-formatted \code{MultiAssayExperiment} object
#' with major OTUs retained.
#' @export
#' @importFrom rlang .data
#' @examples
#'   in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |>
#'     readRDS()
#'   filter_MAE(in_dat, relabu_threshold = 3, occur_pct_cutoff = 5,
#'              taxon_level = "genus")
#'

filter_MAE <- function(dat, relabu_threshold = 3, occur_pct_cutoff = 5,
                       taxon_level = "genus") {
  if(relabu_threshold < 0 | relabu_threshold > 100) {
    stop("relabu_threshold must be between 0 and 100.")
  } else if(occur_pct_cutoff < 0 | occur_pct_cutoff > 100) {
    stop("occur_pct_cutoff must be between 0 and 100.")
  }
  
  # Extract metadata, taxonomic info, and counts
  parsed <- parse_MAE_SE(dat, which_assay = "MicrobeGenetics", type = "MAE")
  tax_table <- parsed$tax
  sam_table <- parsed$sam
  counts_table <- parsed$counts
  # upsample if needed
  up_counts <- animalcules::upsample_counts(counts_table, tax_table,
                                            taxon_level)
  
  samp_range <- range(colSums(counts_table))
  message("The overall range of relative abundance ",
          "counts between samples is (", samp_range[1], ", ",
          samp_range[2], ") \n")
  
  otu.pct <- up_counts %>%
    as.matrix() %>%
    # Get rel abu within samples
    prop.table(margin = 2) 
  
  # prune the OTUs
  # quantile(round(rowMeans(otu.pct), 2))
  otu.select <- apply(otu.pct, 1,
                      # Proportion of samples that the OTU relabu proportion
                      # is greater than relative abundance threshold
                      function(x) sum(x >= relabu_threshold/100) / length(x)) %>%
    # Taking that proportion and checking if that
    # occurrence rate is greater than or equal to the cutoff
    magrittr::is_weakly_greater_than(occur_pct_cutoff / 100)
  # Replace any NAs that were already in table for specified tax level
  tax_table[, taxon_level] <- replace(tax_table[, taxon_level],
                                      is.na(tax_table[, taxon_level]),
                                      "Unknown")
  # Everything to the right of taxon_level should change to "other"
  tax_in_table <- colnames(tax_table)
  tax_ind <- which(tax_in_table == taxon_level)
  tax_to_change <- unique(tax_in_table[tax_ind:length(tax_in_table)])
  tax_table_other <- tax_table |>
    dplyr::mutate(dplyr::across(dplyr::all_of(tax_to_change),
                                function(x) replace(x, !otu.select, 
                                                    "Other")))
  
  # Resum the counts table at lowest level
  low_tax <- tax_in_table[length(tax_in_table)]
  counts_other <- counts_table %>%
    dplyr::bind_cols("OTU" = tax_table_other[, low_tax]) %>%
    dplyr::relocate("OTU") %>%
    dplyr::group_by(.data$OTU) %>%
    dplyr::summarise(dplyr::across(.cols = dplyr::all_of(colnames(counts_table)),
                                   .fns = sum)) %>%
    tibble::column_to_rownames("OTU")
  
  # Adjust tax table accordingly
  tax_table_other2 <- tax_table_other %>%
    dplyr::distinct(.data[[low_tax]], .keep_all = TRUE) %>%
    dplyr::arrange(.data[[low_tax]])
  
  rownames(tax_table_other2) <- tax_table_other2[, low_tax]
  message("Number of OTUs that exhibit a relative abundance >",
          relabu_threshold, "% in at least ",
          occur_pct_cutoff, "% of the total samples: ",
          nrow(tax_table_other2), "/", nrow(tax_table), "\n")
  ## Create SE object
  MAE_out <- create_formatted_MAE(counts_other, tax_table_other2, sam_table)
  return(MAE_out)
}
