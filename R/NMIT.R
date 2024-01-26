# Helper function to rearrange MAE by unit, then time
.rearrange_sub_time <- function(dat, unit_var, covariate_time) {
  all_dat <- parse_MAE_SE(dat)
  map <- all_dat$sam
  otu <- all_dat$counts
  tax <- all_dat$tax
  map2 <- map |>
    dplyr::arrange(!!! rlang::syms(c(unit_var, covariate_time)))
  ind <- match(rownames(map2), colnames(otu))
  otu2 <- otu[, ind]
  MAE <- create_formatted_MAE(otu2, tax, map2)
  return(MAE)
}

#' Calculate within-subject OTU correlations
#'
#' This function takes a \code{MultiAssayExperiment} and outputs an array of
#' temporal intra-subject correlation matrices.
#'
#' @author Yilong Zhang, Huilin Li, Aubrey Odom
#'
#' @inheritParams clean_animalcules_MAE
#' @param method an option of the correlation method ("pearson", "kendall",
#'   "spearman"). The default method is "kendall".
#' @param unit_var a numeric vector of subject.
#' @param fill_na a number between 0 and 1 to fill the missing value. The
#'   default value is 0.
#'
#' @return An three-dimensional array of temporal correlation matrices for each
#'   subject.
#' @export
#' @examples
#' dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
#' output <- tscor(dat, unit_var = "Subject", method = "spearman")
#' head(output)
#' 

tscor <- function(dat, unit_var, method = "kendall", fill_na = 0) {
  if (!(method %in% c("pearson","kendall","spearman"))) {
    stop("method must be one of 'pearson', 'kendall', 'spearman'")
  }
  all_dat <- parse_MAE_SE(dat)
  map <- all_dat$sam
  otus <- all_dat$counts |> t() |> as.data.frame()
  d <- ncol(otus)
  
  # Split OTU tables by unit_var
  sub_split <- map |> dplyr::pull(dplyr::all_of(unit_var))
  otu_list <- split(otus, sub_split)
  
  # Pairwise correlations for each OTU sub-table
  cor_list <- lapply(otu_list, stats::cor, method = method)
  names.cor <- names(cor_list)
  names.otu <- colnames(otus)
  
  # Collect correlations into multidimensional array
  cor_unlist <- unlist(cor_list)
  cor_unlist[is.na(cor_unlist)] <- fill_na
  cor_arr <- array(cor_unlist,
                   dim = c(d, d, length(cor_list)),
                   dimnames = list(names.otu, names.otu, names.cor))
  return(cor_arr)
}

#' Nonparametric Microbial Interdependence Test (NMIT)
#' 
#' The package performs a multivariate distance-based test for group comparisons of
#' microbial temporal interdependence. The NMIT test provides a comprehensive way to evaluate
#' the association between key phenotypic variables and microbial interdependence.
#' This function is recommended for use after a filtering step using \code{filter_animalcules_MAE}.
#'
#' @author Yilong Zhang, Huilin Li, Aubrey Odom
#'
#' @inheritParams tscor
#' @inheritParams plot_stacked_bar
#' @param fixed_cov A character vector of covariates of interest found in dat.
#' @param dist_type A character strin specifying the type of matrix norm to be computed. The default is \code{"F"}. 
#' \itemize{
#'     \item{\code{"M"} or \code{"m"}}{specifies the maximum modulus of all the elements in \code{x}.}
#'     \item{\code{"O"}, \code{"o"} or \code{"1"}}{specifies the one norm, (maximum absolute column sum);}
#'     \item{\code{"I"} or \code{"i"}}{specifies the infinity norm (maximum absolute row sum);}
#'     \item{\code{"F"} or \code{"f"}}{specifies the Frobenius norm (the Euclidean norm of \code{x} treated as if it were a vector); and}
#' }  
#' @param heatmap A logical value indicating whether to draw heatmap. The default is \code{TRUE}.
#' @param classify A logical value indicating whether to draw a classifier tree. The default is \code{FALSE}.
#' @param fill_na A number between 0 and 1 to fill \code{NA} values. The default value is 0.
#'
#' @return This function returns an analysis of variance (ANOVA) table showing sources of variation,
#' degrees of freedom, sequential sums of squares, mean squares, F statistics, partial R-squared and P values,
#' based on 999 permutations.
#' 
#' @export
#' 
#' @examples
#' dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
#' NMIT(dat, unit_var = "Subject", fixed_cov = "Group", covariate_time = "Month")
#'

NMIT <- function(dat, unit_var, fixed_cov,
                 covariate_time,
                 method = "kendall", dist_type = "F",
                 heatmap = TRUE, classify = FALSE,
                 fill_na = 0) {
  dat <- .rearrange_sub_time(dat, unit_var, covariate_time) 
  all_dat <- parse_MAE_SE(dat)
  map <- all_dat$sam
  otu <- all_dat$counts
  map.unique <- unique(map[, c(unit_var, fixed_cov)])
  rownames(map.unique) <- map.unique[, unit_var]
  n.sample <- length(unique(map[, unit_var]))
  otu.cor <- tscor(dat, method = method,
                   unit_var = unit_var, fill_na = fill_na) |>
    suppressWarnings()
  dist <- outer(seq_len(n.sample), seq_len(n.sample),
                function(x,y) {
                  foo <- function(x,y) norm(otu.cor[, , x] - otu.cor[, , y],
                                            type = dist_type)
                  mapply(foo, x,y)
                })
  rownames(dist) <- dimnames(otu.cor)[[3]]
  colnames(dist) <- dimnames(otu.cor)[[3]]
  grp <- map.unique[rownames(dist) , fixed_cov]
  test <- vegan::adonis2(dist ~ grp)
  if(heatmap){
    #pvalue <- paste(round(c(test$F[1], test$`Pr(>F)`[1]), 3),
    #                collapse = ";")
    pvalue <- round(c(test$`Pr(>F)`[1]), 3)
    n.taxa <- nrow(otu)
    diag(dist) <- stats::median(dist)
    plot_title <- paste0(" #OTU:", n.taxa, " pvalue:", pvalue)
    stats::heatmap(sqrt(dist), Colv = classify, Rowv = classify, scale = "column", 
            xlab = unit_var, ylab = unit_var,
            main = paste("NMIT Correlation Heatmap across", n.taxa, "taxa with p-value", pvalue))
  }
}
