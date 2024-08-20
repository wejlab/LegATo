# Helper function: unpaired multivariate Hotelling's t-squared

.Hotelling_mv_T2_un <- function(input_data, Populations, Subjects,
                                taxon = "taxon", save_table_loc, num_taxa) {
  # Rename groups
  input_data <- input_data %>%
    dplyr::select("Populations" = dplyr::starts_with(Populations),
                  "Subjects" = dplyr::starts_with(Subjects),
                  "Taxon" = dplyr::starts_with(taxon),
                  "Abundance")
  # Check that there is only one observation per group.
  nonuniq_obs <- input_data %>%
    dplyr::group_by(.data$Populations, .data$Subjects, .data$Taxon) %>%
    dplyr::summarize("count" = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(.data$count > 1) %>%
    nrow()
  if(nonuniq_obs > 0) message("More than one observation per unit/group detected") 
  # Back to groups
  Group1 <- unique(input_data$Populations)[1]
  Group2 <- unique(input_data$Populations)[2]
  Sub1 <- input_data %>% dplyr::filter(.data$Populations == Group1) %>% 
    dplyr::distinct(.data$Subjects) %>% dplyr::ungroup() %>% 
    dplyr::select("Subjects") %>% unlist()
  Sub2 <- input_data %>% dplyr::filter(.data$Populations == Group2) %>% 
    dplyr::distinct(.data$Subjects) %>% dplyr::ungroup() %>%
    dplyr::select("Subjects") %>% unlist()
  # Define n
  n <- input_data %>% dplyr::group_by(Populations) %>%
    dplyr::distinct(Subjects) %>% dplyr::summarize("n_col" = dplyr::n())
  n1 <- n %>% dplyr::filter(Populations == Group1) %>%
    dplyr::select("n_col") %>% as.numeric()
  n2 <- n %>% dplyr::filter(Populations == Group2) %>%
    dplyr::select("n_col") %>% as.numeric()
  p <- length(unique(input_data$Taxon))
  # Check that it won't be singular
  # (n1 + n2 - num_taxa - 1) < 1
  num_taxa_min <- n1 + n2 - 2
  if (num_taxa > num_taxa_min) stop("num_taxa must be no larger than ", num_taxa_min)
  # Sample mean vector
  X_i <- input_data %>%
    dplyr::rename("Xi" = "Abundance")
  Xbar <- X_i %>%
    dplyr::group_by(.data$Populations, .data$Taxon) %>%
    dplyr::reframe(Xbar = mean(.data$Xi))
  # X_i - Xbar
  diff <- X_i %>%
    dplyr::left_join(., Xbar, by = c("Populations", "Taxon")) %>%
    dplyr::group_by(.data$Subjects, .data$Taxon) %>%
    dplyr::reframe("diff" = .data$Xi - .data$Xbar)
  # (X_ij-Xbar_i)%*%(X_ij-Xbar_i)'
  mult_func <- function(x) {
    vec <- t(t(diff$diff[diff$Subjects == x]))
    vec %*% t(vec)
  }
  # Calculate S_p
  # Sample Var-cov matrix of vecs Y_i
  S_1 <- 1/(n1 - 1) * base::Reduce("+", lapply(Sub1, mult_func))
  S_2 <- 1/(n2 - 1) * base::Reduce("+", lapply(Sub2, mult_func))
  S_p <- ((n1 - 1) * S_1 + (n2 - 1) * S_2) / (n1 + n2 - 2)
  
  # T^2 = t(Xbar_1-Xbar_2) %*% {S_p(1/n1 + 1/n2)}^-1 %*% (Xbar_1-Xbar_2)
  Xbar_1 <- Xbar %>% dplyr::filter(.data$Populations == Group1) %>% dplyr::select("Xbar")
  Xbar_2 <- Xbar %>% dplyr::filter(.data$Populations == Group2) %>% dplyr::select("Xbar")
  meandiff <- t(t(Xbar_1 - Xbar_2))
  T_2 <- t(meandiff) %*% solve(S_p * (1 / n1 + 1 / n2)) %*% meandiff
  
  #F statistic
  # Dist according to F_p, n1 + n2 - p - 1
  F_stat <- (n1 + n2 - p - 1) / (p * (n1 + n2 - 2)) * T_2
  # Reject H0 at alpha if F-val exceeds
  df1 <- p
  df2 <- n1 + n2 - p - 1
  crit_F <- stats::qf(0.95, df1, df2)
  pval <- 1 - stats::pf(c(F_stat), df1, df2)
  
  if(pval < 0.05) {
    ttest_unpaired <- function(tax) {
      for_testing <- input_data %>% dplyr::filter(.data$Taxon == tax)
      out_test <- stats::t.test(stats::formula("Abundance ~ Populations"), for_testing,
                                alternative = "two.sided",
                                var.equal = FALSE, paired = FALSE)
      return(list(t = out_test$statistic, df = out_test$parameter,
                  diff_means = out_test$estimate[1],
                  CI_2.5 = out_test$conf.int[1],
                  CI_97.5 = out_test$conf.int[2],
                  "p-value" = out_test$p.value))
    }
    results <- plyr::aaply(unique(input_data$Taxon), 1, ttest_unpaired)
    results <- rbind(results,
                     "adj p-value" = stats::p.adjust(results["p-value", ],
                                                     method = "bonferroni"))
    # Write results
    filename_out <- paste0("ttest_results", Group1, Group2, ".csv")
    utils::write.csv(results, file.path(save_table_loc, filename_out))
    message("Results from t-tests written to", filename_out)
  }
  
  list("df1" = df1, "df2" = df2, "crit_F" = crit_F,
       "F_stat" = F_stat[1], "pvalue" = pval)
}

# Helper function: Paired Multivariate Hotelling's T-Squared Test

.Hotelling_mv_T2 <- function(input, Group1, Group2, save_table_loc){
  # Rename groups
  input_data <- input %>% dplyr::select("pairing",
                                        "Taxon" = "taxon",
                                        "Group1" = tidyr::starts_with(Group1),
                                        "Group2" = tidyr::starts_with(Group2))
  # Check that there is only one observation per group.
  nonuniq_obs <- input_data %>%
    dplyr::group_by(.data$pairing, .data$Taxon) %>%
    dplyr::summarize("count" = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(.data$count > 1) %>%
    nrow()
  if(nonuniq_obs > 0) message("More than one observation per unit/group detected")
  # Define n
  all_pairs <- unique(input_data$pairing)
  n <- length(all_pairs)
  p <- length(unique(input_data$Taxon))
  # Check that it won't be singular
  # n - p < 1
  if (n - p < 1) stop("num_taxa must be no larger than ", n - 1)
  # Sample mean vector
  Y_i <- input_data %>%
    dplyr::group_by(.data$pairing, .data$Taxon) %>%
    dplyr::summarise(Yi = .data$Group1 - .data$Group2,
                     .groups = "drop")
  Ybar <- Y_i %>%
    dplyr::group_by(.data$Taxon) %>%
    dplyr::summarise(Ybar = mean(.data$Yi),
                     .groups = "drop")
  # Y_i - Ybar
  diff <- Y_i %>%
    dplyr::left_join(., Ybar, by = "Taxon") %>%
    dplyr::group_by(.data$pairing, .data$Taxon) %>%
    dplyr::summarise("diff" = .data$Yi - .data$Ybar, .groups = "drop")
  # (y_i-Ybar)%*%(y_i-Ybar)'
  mult_func <- function(x) {
    vec <- t(t(diff$diff[diff$pairing == x]))
    vec %*% t(vec)
  }
  # Calculate S_Y
  # Sample Var-cov matrix of vecs Y_i
  S_Y <- 1/(n - 1) * base::Reduce("+", lapply(all_pairs, mult_func))
  # Calculate T^2
  # T^2 = n*Ybar'*inv(S_Y)*Ybar
  Ybar_v <- t(t(Ybar$Ybar))
  T_2 <- n * t(Ybar_v) %*% solve(S_Y) %*% Ybar_v
  #F statistic
  # Dist according to F_p,n-p
  F_stat <- (n - p) / (p * (n - 1)) * T_2
  # Reject H0 at alpha if F-val exceeds
  df1 <- p
  df2 <- n - p
  crit_F <- stats::qf(0.95, df1, df2)
  pval <- 1 - stats::pf(c(F_stat), df1, df2)
  
  # Conduct a paired t-test
  # Following # 1 in R&C pg. 140
  if(pval < 0.05) {
    ttest_paired <- function(tax) {
      for_testing <- input_data %>% dplyr::filter(.data$Taxon == tax)
      out_test <- stats::t.test(x = for_testing$Group1,
                                y = for_testing$Group2,
                                alternative = "two.sided",
                                paired = TRUE, var.equal = TRUE)
      output <- list("t" = as.numeric(out_test$statistic),
                     "df" = as.numeric(out_test$parameter),
                     "diff_means" = as.numeric(out_test$estimate),
                     "CI_2.5" = as.numeric(out_test$conf.int[1]),
                     "CI_97.5" = as.numeric(out_test$conf.int[2]),
                     "p-value" = as.numeric(out_test$p.value))
      return(lapply(output, base::round, digits = 4))
    }
    results <- plyr::aaply(unique(input_data$Taxon), 1, ttest_paired)
    results_2 <- cbind(results, "adj p-value" = stats::p.adjust(results[, "p-value"],
                                                                method = "bonferroni"))
    # Write results
    filename_out <- paste0("ttest_results", Group1, Group2, ".csv")
    utils::write.csv(results_2, file.path(save_table_loc, filename_out))
    message("Results from t-tests written to ", filename_out)
  }
  
  return(list("df1" = df1, "df2" = df2, "crit_F" = crit_F,
              "F_stat" = F_stat[1], "pvalue" = pval))
}

#' Conduct a multivariate Hotelling's T-squared test
#'
#' This function takes an animalcules-formatted \code{MultiAssayExperiment}
#' object and runs a multivariate Hotelling's T-squared test. The test expects
#' a comparison of two distinct groups, and compares the abundances of the top microbes at a
#' given taxon level between the groups. This function allows both paired and unpaired
#' tests. Both test the null hypothesis that the population
#' mean vectors are equal, with the alternative being that they are unequal.
#' 
#' The Hotelling's t-squared statistic (t2) is a generalization of
#' Student's t-statistic that is used in multivariate hypothesis testing
#' to test the means of different populations.
#' 
#' Note that any entries or pairs with missing values are excluded.
#' 
#' Referenced articles in the implementation of tests:
#'
#' https://online.stat.psu.edu/stat505/lesson/7/7.1/7.1.14
#' 
#' https://online.stat.psu.edu/stat505/lesson/7/7.1/7.1.15
#' 
#' https://online.stat.psu.edu/stat505/lesson/7/7.1/7.1.4
#' 
#' https://online.stat.psu.edu/stat505/lesson/7/7.1/7.1.9
#'
#' @param dat A MultiAssayExperiment object specially formatted as an animalcules output.
#' @param test_index Any argument used for subsetting the input \code{dat},
#' can be a character, logical, integer, list or List vector. Default is \code{NULL}.
#' @param taxon_level Character string, default is \code{"genus"}.
#' @param num_taxa The number of most abundant taxa to test. If unpaired, this should
#' be no larger than the total number of subjects in both groups - 2, or (n1 + n2 -2).
#' If paired, this should be no larger than the total number of pairs - 1, or n - 1.
#' Required.
#' @param grouping_var Character string, the name of a DICHOTOMOUS grouping variable in the
#' metadata of \code{dat}.
#' @param paired Logical indicating whether a paired test should be conducted.
#' Default is \code{FALSE} for an unpaired test.
#' @param pairing_var Character string giving the variable containing pairing information.
#' The variable should be in integer form. Must be supplied if \code{paired = TRUE},
#' otherwise the default is NULL.
#' @param unit_var Character string giving the variable containing the identifiers
#' for the unit on which multiple measurements were conducted, e.g. subjects. Default is
#' \code{NULL}; must be supplied if \code{paired = FALSE}.
#' @param save_table_loc A character string giving the folder path to save t.test results.
#' Note that these are only conducted if the Hotelling's T-test value is <0.05. 
#' Defaults to the current working directory.
#' 
#' @return A list of the elements "df1", "df2", "crit_F", "F_stat" and "pvalue" giving the
#' results of the test.
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' dat <- system.file("extdata", "MAE.RDS", package = "LegATo") |>
#' readRDS()
#' dat_0.05 <- filter_MAE(dat)
#' out1 <- test_hotelling_t2(dat = dat_0.05,
#'                   test_index = which(dat_0.05$MothChild == "Infant" &
#'                                        dat_0.05$timepoint == 6),
#'                   taxon_level = "genus",
#'                   num_taxa = 6,
#'                   paired = TRUE,
#'                   grouping_var = "HIVStatus",
#'                   pairing_var = "pairing")
#' out1                  
#' 
#' out <- test_hotelling_t2(dat = dat_0.05,
#'                   test_index = which(dat_0.05$MothChild == "Mother" &
#'                                        dat_0.05$timepoint == 0),
#'                   taxon_level = "genus",
#'                   num_taxa = 12,
#'                   grouping_var = "HIVStatus",
#'                   unit_var = "Subject",
#'                   paired = FALSE)
#' out                  
#' 

test_hotelling_t2 <- function(dat,
                              test_index = NULL,
                              taxon_level = "genus",
                              num_taxa, grouping_var,
                              paired = FALSE,
                              pairing_var = NULL,
                              unit_var = NULL,
                              save_table_loc = ".") {
  MAE <- dat
  if (!is.null(test_index)) MAE <- MultiAssayExperiment::subsetByColData(dat, test_index)
  best_taxa <- get_top_taxa(MAE, taxon_level) %>%
    dplyr::pull(.data$taxon)
  input_data <- get_long_data(MAE, taxon_level, log = TRUE,
                              counts_to_CPM = TRUE) %>%
    dplyr::filter(.data$taxon %in% best_taxa[seq_len(num_taxa)]) 
  
  if (paired) {
    if (is.null(pairing_var)) message("Please supply pairing_var if paired = TRUE")
    output_data <- input_data %>%
      dplyr::select(dplyr::all_of(c(grouping_var, pairing_var)), "taxon",
                    "Abundance") %>%
      tidyr::pivot_wider(., id_cols = dplyr::all_of(c("taxon", pairing_var)),
                         values_from = "Abundance",
                         names_from = dplyr::all_of(grouping_var)) %>%
      dplyr::arrange(.data$pairing) %>%
      dplyr::filter(!dplyr::if_any(tidyselect::everything(), is.na))
    # Ensure dichotomous
    all_vars <- unique(input_data[, grouping_var])
    if(length(all_vars) != 2) {
      stop("grouping_var must be dichotomous")
    } else results <- .Hotelling_mv_T2(output_data, all_vars[1],
                                       all_vars[2], save_table_loc)
  } else {
    if (is.null(unit_var)) message("Please supply unit_var if paired = FALSE")
    output_data <- input_data %>%
      dplyr::select("taxon", dplyr::all_of(c(unit_var, grouping_var)), "Abundance") %>%
      dplyr::arrange(.data[[unit_var]]) %>%
      dplyr::filter(!dplyr::if_any(tidyselect::everything(), is.na))
    results <- .Hotelling_mv_T2_un(output_data, grouping_var, unit_var, "taxon",
                                   save_table_loc, num_taxa)
  }
  return(results)
}
