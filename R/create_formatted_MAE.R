#' Create a formatted MultiAssayExperiment compatible with LegATo
#'
#' This function takes either a counts_dat, tax_dat, and metadata_dat input OR
#' a TreeSummarizedExperiment input and creates
#' a specifically-formatted MAE object that is compatible for use with LegATo
#' and animalcules. Checks are performed on inputs to ensure that they can be
#' integrated properly.
#'
#' @param counts_dat A matrix, data.table, or data.frame consisting of microbial
#'   raw counts data. The \code{colnames} should be sample names and the
#'   \code{rownames} should be in the same order as the \code{tax_dat} entries.
#'   Not required if \code{tree_SE} is passed in.
#' @param tax_dat A matrix, data.table, or data.frame of hierarchical taxonomic
#'   data. Should have columns such as "family", "genus", "species" with each
#'   row uniquely delineating a different taxon. The rows should be in the same
#'   order as the rows of \code{counts_dat}. Not required if \code{tree_SE} is
#'   passed in.
#' @param metadata_dat A metadata table with \code{rownames} equivalent to the
#'   samples that are the \code{colnames} of the \code{counts_dat}. Not required
#'   if \code{tree_SE} is passed in.
#' @param tree_SE A TreeSummarizedExperiment object with counts, taxonomy, and
#' metadata.
#'   
#' @export
#' @returns A \code{MultiAssayExperiment} object.
#' @examples
#' nsample <- ntaxa <- 3
#' counts_dat <- data.frame(
#'     "X123" = runif(ntaxa, 0, 500),
#'     "X456" = runif(ntaxa, 0, 500),
#'     "X789" = runif(ntaxa, 0, 500)
#' )
#' tax_dat <- data.frame(
#'     "class" = c("rand1", "rand2", "rand3"),
#'     "species" = c("rand4", "rand5", "rand6")
#' ) |>
#'     as.data.frame()
#' # Set rownames as lowest unique taxonomic level
#' rownames(tax_dat) <- tax_dat$species
#' rownames(counts_dat) <- tax_dat$species
#' metadata <- data.frame(
#'     Sample = c("X123", "X456", "X789"),
#'     Group = c("A", "B", "A"),
#'     Var = rnorm(nsample)
#' )
#' rownames(metadata) <- metadata$Sample
#' out_MAE <- create_formatted_MAE(counts_dat, tax_dat, metadata)
#' 
#' # TreeSummarizedExperiment
#' tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
#'     assays = list(counts = counts_dat),
#'     colData = metadata,            
#'     rowData = tax_dat             
#' )
#' out_MAE_2 <- create_formatted_MAE(tree_SE = tse)
#' out_MAE_2

create_formatted_MAE <- function(counts_dat = NULL, tax_dat = NULL, metadata_dat = NULL,
                                 tree_SE = NULL) {
    if (methods::is(tree_SE, "TreeSummarizedExperiment")) {
        counts_dat <- SummarizedExperiment::assays(tree_SE)[[1]]
        tax_dat <-SummarizedExperiment::rowData(tree_SE)
        metadata_dat <- SummarizedExperiment::colData(tree_SE) |> as.data.frame() 
    } else {
        if (is.null(counts_dat) | is.null(tax_dat) | is.null(metadata_dat)) {
            stop("Please supply counts, taxonomy, and metadata tables.")
        }
        to_check <- c("matrix", "data.table", "data.frame", "tibble")
        if (!(class(counts_dat) %in% to_check)) stop("counts_dat should be one of",
                                                 to_check)
        if (!(class(metadata_dat) %in% to_check)) stop("metadata_dat should be one of",
                                                   to_check)
        if (!(class(tax_dat) %in% to_check)) stop("tax_dat should be one of",
                                              to_check)
        if(nrow(tax_dat) != nrow(counts_dat)) stop(
            "The number of rows of tax_dat and counts_dat should be equal.")
        if(nrow(metadata_dat) != ncol(counts_dat)) stop(
            "The number of rows of metadata_dat and columns of counts_dat should be equal.")
    }
    
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
