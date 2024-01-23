get_top_taxa <- function(dat, taxa = "genus") {
  microbe <- parse_MAE_SE(dat)
  all_relabu <- microbe$counts |>
    animalcules::upsample_counts(microbe$tax, taxa)%>%
    animalcules::counts_to_relabu() %>%
    tibble::rownames_to_column(var = "taxa") |>
    dplyr::rowwise(taxa) |>
    # Sum everything but the first columm ("genus")
    dplyr::summarise(
      allmeans = mean(
        dplyr::c_across(dplyr::all_of(colnames(microbe$counts)))),
      .groups = "drop") |>
    dplyr::arrange(desc(allmeans))
  
  return(all_relabu)
}
