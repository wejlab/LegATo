utils::globalVariables(".")

#' Filter a MultiAssayExperiment to a top percentage of taxa and label the rest
#' as "Other"
#'
#' This function takes an animalcules-formatted \code{MultiAssayExperiment}
#' (MAE) object and identifies all taxa at the "genus" level that represent
#' <\code{filter_prop} average relative abundance across all samples in the MAE.
#' After identification at the genus level, taxa across the genus and species
#' levels are then consolidated into the category "Other".
#'
#' @inheritParams plot_stacked_bar
#' @param filter_prop A double strictly between 0 and 1, representing
#' the proportion of relative abundance at which to filter. Default is 0.001.
#'
#' @returns An animalcules-formatted \code{MultiAssayExperiment} object with
#'   appropriate filtration.
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
#' filter_animalcules_MAE(in_dat, 0.01)
#' 

filter_animalcules_MAE <- function(dat, filter_prop = 0.001) {
    if(filter_prop <= 0 | filter_prop >= 1) stop("filter_prop must be between 0 and 1.")
    # Extract metadata, taxonomic info, and counts
    parsed <- parse_MAE_SE(dat, which_assay = "MicrobeGenetics", type = "MAE")
    tax_table <- parsed$tax
    sam_table <- parsed$sam
    counts_table <- parsed$counts
    
    all_relabu_genus <- counts_table %>%
        as.matrix() %>%
        # Get rel abu within samples
        prop.table(margin = 2) %>%
        tibble::as_tibble() %>%
        dplyr::bind_cols("genus" = tax_table$genus) %>%
        dplyr::relocate("genus") %>%
        dplyr::group_by(.data$genus) %>%
        # Sum rel abu within samples/columns for genera
        dplyr::summarise(dplyr::across(.fns = sum, .cols = dplyr::everything())) %>%
        # Sum everything but the first columm ("genus")
        dplyr::mutate("allmeans" = rowMeans(dplyr::select(., -1))) %>%
        dplyr::select("genus", "allmeans") %>%
        dplyr::mutate("genus" = replace(.data$genus, is.na(.data$genus), "Unknown")) %>%
        dplyr::arrange(dplyr::desc(.data$allmeans)) %>%
        dplyr::mutate("lessthanXpct" = .data$allmeans < filter_prop)
    
    # Identify species in other
    othergenera <- all_relabu_genus %>%
        dplyr::filter(.data$lessthanXpct == TRUE) %>%
        dplyr::select("genus") %>% unlist() %>% unname()
    
    ## Use the identified species to update the tax, species tables
    
    # Replace tax table
    tax_table_other <- tax_table %>%
        dplyr::mutate("genus" = replace(.data$genus, is.na(.data$genus), "Unknown")) %>%
        dplyr::mutate(dplyr::across(c(.data$genus, .data$species),
                                    function(x) replace(x, .data$genus %in% othergenera, "Other")))
    
    # Resum the counts table
    counts_other <- counts_table %>%
        dplyr::bind_cols("species" = tax_table_other$species) %>%
        dplyr::group_by(.data$species) %>%
        dplyr::summarise(dplyr::across(.cols = dplyr::all_of(colnames(counts_table)),
                                       .fns = sum)) %>%
        tibble::column_to_rownames("species")
    
    # Adjust tax table accordingly
    tax_table_other2 <- tax_table_other %>%
        dplyr::distinct(.data$species, .keep_all = TRUE) %>%
        dplyr::arrange(.data$species)
    
    rownames(tax_table_other2) <- tax_table_other2$species
    
    ## Create SE object
    MAE_out <- create_formatted_MAE(counts_other, tax_table_other2, sam_table)
    return(MAE_out)
}
