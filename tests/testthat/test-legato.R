

testthat::test_that("plot_alluvial works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    testthat::expect_s3_class(
        plot_alluvial(
            in_dat,
            taxon_level = "family",
            covariate_1 = "Group",
            covariate_time = "Month",
            palette_input = rainbow(25)
        ),
        "ggplot"
    )
    
})

testthat::test_that("plot_spaghetti works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    all_taxa <- get_top_taxa(in_dat, "phylum")
    
    testthat::expect_s3_class(
        plot_spaghetti(
            in_dat,
            taxon_level = "phylum",
            covariate_1 = "Group",
            covariate_time = "Month",
            unit_var = "Subject",
            which_taxon = all_taxa$taxon[1],
            palette_input = rainbow(25)
        ),
        "ggplot"
    )
    
})

testthat::test_that("plot_stacked_area works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    
    testthat::expect_s3_class(
        plot_stacked_area(
            in_dat,
            taxon_level = "phylum",
            covariate_1 = "Group",
            covariate_time = "Month",
            palette_input = rainbow(25)
        ),
        "ggplot"
    )
})

testthat::test_that("plot_stacked_bar works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    
    testthat::expect_s3_class(
        plot_stacked_bar(
            in_dat,
            taxon_level = "family",
            covariate_1 = "Group",
            covariate_time = "Month",
            palette_input = rainbow(25)
        ),
        "ggplot"
    )
})

testthat::test_that("plot_heatmap works", {
    mat_testdata <- rbind(matrix(c(rnorm(80), rnorm(80) + 5), 16, 10, dimnames = list(
        paste0("Taxon", seq_len(16)), paste0("sample", seq_len(10))
    )),
    matrix(rnorm(1000), 100, 10, dimnames = list(
        paste0("Taxon0", seq_len(100)), paste0("sample", seq_len(10))
    )))
    cov_mat <- data.frame(sample = c(rep("down", 5), rep("up", 5))) |>
        magrittr::set_rownames(paste0("sample", seq_len(10)))
    
    # Example using custom colors for the annotation information
    color2 <- stats::setNames(c("purple", "black"), c("down", "up"))
    color.list <- list("sample" = color2)
    
    testthat::expect_type(
        plot_heatmap(
            inputData = mat_testdata,
            annotationData = cov_mat,
            name = "Data",
            plot_title = "Example",
            plottingColNames = NULL,
            annotationColNames = NULL,
            colList = color.list,
            scale = FALSE,
            showColumnNames = TRUE,
            showRowNames = FALSE,
            colorSets = c(
                "Set1",
                "Set2",
                "Set3",
                "Pastel1",
                "Pastel2",
                "Accent",
                "Dark2",
                "Paired"
            ),
            choose_color = c("blue", "gray95", "red"),
            split_heatmap = "none",
            column_order = NULL
        ),
        "NULL"
    )
})

testthat::test_that("distinctColors works", {
    testthat::expect_vector(distinctColors(10))
})

testthat::test_that("run_gee_model works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |>
        readRDS()
    testthat::expect_s3_class(
        run_gee_model(
            in_dat,
            taxon_level = "genus",
            unit_var = "Subject",
            fixed_cov = c("HairLength", "Age", "Group", "Sex"),
            corstr = "ar1"
        ),
        "data.frame"
    )
})

testthat::test_that("filter_animalcules_MAE works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    
    
    testthat::expect_s4_class(filter_animalcules_MAE(in_dat, 0.01),
                              "MultiAssayExperiment")
})

testthat::test_that("run_lmm_model works", {
   dat <- system.file("extdata/MAE.RDS", package = "LegATo") |>
       readRDS() |>
       filter_animalcules_MAE(0.05)
   testthat::expect_s3_class(run_lmm_model(
       dat,
       taxon_level = "genus",
       unit_var = "Subject",
       fixed_cov = c("HIVStatus", "timepoint")
   ), "data.frame")
})

testthat::test_that(
    "run_lm_model works",
    {
        dat <- system.file("extdata/MAE.RDS", package = "LegATo") |>
            readRDS() |>
            filter_MAE(0.05)
        testthat::expect_s3_class(run_lm_model(
            dat,
            fixed_cov = c("timepoint", "HIVStatus"),
            plot_out = FALSE
        ),
        "data.frame")
})

testthat::test_that("test_hotelling_t2 works", {
    dat <- system.file("extdata", "MAE.RDS", package = "LegATo") |>
        readRDS()
    dat_0.05 <- filter_MAE(dat, 0.05)
    
    testthat::expect_type(
        test_hotelling_t2(
            dat = dat_0.05,
            test_index = which(dat_0.05$MothChild == "Infant" &
                                   dat_0.05$timepoint == 6),
            taxon_level = "genus",
            # To avoid n < p, use top 5-6 species
            num_taxa = 6,
            paired = TRUE,
            grouping_var = "HIVStatus",
            pairing_var = "pairing"
        ),
        "list"
    )
    testthat::expect_type(
        test_hotelling_t2(
            dat = dat_0.05,
            test_index = which(dat_0.05$MothChild == "Mother" &
                                   dat_0.05$timepoint == 6),
            taxon_level = "genus",
            # To avoid n < p, use top 5-6 species
            num_taxa = 6,
            paired = FALSE,
            unit_var = "Subject",
            grouping_var = "HIVStatus"
        ),
        "list"
    )
})

testthat::test_that("NMIT works", {
    dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    testthat::expect_type(NMIT(
        dat,
        unit_var = "Subject",
        fixed_cov = "Group",
        covariate_time = "Month"
    ), "NULL")
})

 testthat::test_that("tscor works", {
     dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |>
         readRDS() |>
         filter_animalcules_MAE(0.05)
     testthat::expect_type(tscor(dat, unit_var = "Subject",
                                 method = "spearman", fill_na = 0.0004), "double")
 })

testthat::test_that("create_formatted_MAE works", {
    nsample <- ntaxa <- 3
    counts_dat <- data.frame(
        "X123" = runif(ntaxa, 0, 500),
        "X456" = runif(ntaxa, 0, 500),
        "X789" = runif(ntaxa, 0, 500)
    )
    tax_dat <- data.frame(
        "class" = c("rand1", "rand2", "rand3"),
        "species" = c("rand4", "rand5", "rand6")
    ) |>
        as.data.frame()
    # Set rownames as lowest unique taxonomic level
    rownames(tax_dat) <- tax_dat$species
    rownames(counts_dat) <- tax_dat$species
    metadata <- data.frame(
        Sample = c("X123", "X456", "X789"),
        Group = c("A", "B", "A"),
        Var = rnorm(nsample)
    )
    rownames(metadata) <- metadata$Sample
    testthat::expect_s4_class(create_formatted_MAE(counts_dat, tax_dat, metadata),
                              "MultiAssayExperiment")
    # TreeSummarizedExperiment
    tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = list(counts = counts_dat),
        colData = metadata,            
        rowData = tax_dat             
    )
    testthat::expect_s4_class(create_formatted_MAE(tree_SE = tse),
                              "MultiAssayExperiment")
})

testthat::test_that("clean_MAE works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    
    
    testthat::expect_s4_class(clean_MAE(in_dat), "MultiAssayExperiment")
})

testthat::test_that("filter_MAE works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    
    testthat::expect_s4_class(filter_MAE(in_dat, 0.01), "MultiAssayExperiment")
})

testthat::test_that("parse_MAE_SE works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    testthat::expect_type(parse_MAE_SE(in_dat), "list")
})

testthat::test_that("get_long_data works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    testthat::expect_s3_class(get_long_data(in_dat, "genus", log = TRUE, counts_to_CPM = TRUE),
                              "data.frame")
})

testthat::test_that("get_stacked_data works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    
    testthat::expect_s3_class(
        get_stacked_data(
            in_dat,
            "genus",
            covariate_1 = "Sex",
            covariate_time = "Month"
        ),
        "data.frame"
    )
    testthat::expect_s3_class(
        get_stacked_data(
            in_dat,
            "genus",
            covariate_1 = NULL,
            covariate_time = "Month"
        ),
        "data.frame"
    )
})

testthat::test_that("get_summary_table works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    testthat::expect_s3_class(get_summary_table(in_dat, c("Group", "Subject")), "data.frame")
})

testthat::test_that("get_top_taxa works", {
    in_dat <- system.file("extdata/MAE_small.RDS", package = "LegATo") |> readRDS()
    testthat::expect_s3_class(get_top_taxa(in_dat, "genus"), "data.frame")
})
