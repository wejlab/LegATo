#' Create a formatted MultiAssayExperiment compatible with LegATo
#'
#' This function takes a counts_dat, tax_dat, or metadata_dat input and creates
#' a specifically-formatted MAE object that is compatible for use with LegATo
#' and animalcules. Checks are performed on inputs to ensure that they can be
#' integrated properly.
#'
#' @param counts_dat A matrix, data.table, or data.frame consisting of microbial
#'   raw counts data. The \code{colnames} should be sample names and the
#'   \code{rownames} should be in the same order as the \code{tax_dat} entries.
#' @param tax_dat A matrix, data.table, or data.frame of hierarchical taxonomic
#'   data. Should have columns such as "family", "genus", "species" with each
#'   row uniquely delineating a different taxon. The rows should be in the same
#'   order as the rows of \code{counts_dat}.
#' @param metadata_dat A metadata table with \code{rownames} equivalent to the
#'   samples that are the \code{colnames} of the \code{counts_dat}.
#'   
#' @export
#' @returns A \code{MultiAssayExperiment} object.
#' @examples
#' nsample <- ntaxa <- 3
#' counts_dat <- data.frame("X123" = runif(ntaxa, 0, 500),
#'                          "X456" = runif(ntaxa, 0, 500),
#'                          "X789" = runif(ntaxa, 0, 500))
#' tax_dat <- data.frame("class" = c("rand1", "rand2", "rand3"),
#'                       "species" = c("rand4", "rand5", "rand6")) |>
#'   as.data.frame()
#' # Set rownames as lowest unique taxonomic level
#' rownames(tax_dat) <- tax_dat$species
#' metadata <- data.frame(Sample = c("X123", "X456", "X789"),
#'                        Group = c("A", "B", "A"),
#'                        Var = rnorm(nsample))
#' rownames(metadata) <- metadata$Sample
#' out_MAE <- create_formatted_MAE(counts_dat, tax_dat, metadata)
#' out_MAE
#'

create_formatted_MAE <- function(counts_dat, tax_dat, metadata_dat) {
  # Check that inputs are matrix, data.table, or data.frame?
  # Check that the row and column names conform, same dimensions...
  # All taxon names should be lowercase
  se_mgx <- counts_dat %>% base::data.matrix() %>% S4Vectors::SimpleList() %>%
    magrittr::set_names("MGX")
  se_rowData <- tax_dat %>% base::data.frame() %>%
    dplyr::mutate_all(as.character) %>% S4Vectors::DataFrame()
  se_colData  <- metadata_dat %>% S4Vectors::DataFrame()
  
  microbe_se <- SummarizedExperiment::SummarizedExperiment(
    assays = se_mgx, colData = se_colData, rowData = se_rowData)
  MAE_out <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList(MicrobeGenetics = microbe_se),
    colData = se_colData)
  return(MAE_out)
}
