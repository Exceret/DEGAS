#' @title Label Binary Phenotype Cells Based on Prediction Scores with Minimum Threshold
#'
#' @description
#' Classifies cells into binary phenotype groups ("Positive" vs "Other") based on
#' prediction score differences between two phenotypic conditions. This enhanced
#' version includes a minimum threshold constraint to ensure biological relevance
#' and provides detailed reporting on classification outcomes.
#'
#' @param pred_dt A data.table containing prediction scores for two phenotypic
#'   conditions. Must contain columns specified in `pheno_colnames`.
#' @param pheno_colnames Character vector of length 2 specifying the column names
#'   for the two phenotypic conditions to compare. The second element is used as
#'   the reference group if not found with regex matching.
#' @param select_fraction Numeric value between 0 and 1 specifying the fraction
#'   of cells to classify as "Positive". Default selection depends on the
#'   distribution characteristics.
#' @param test_method Character string specifying the statistical test to use
#'   for normality assessment. One of: `"jarque-bera"`, `"d'agostino"`,
#'   `"kolmogorov-smirnov"`.
#' @param min_threshold Numeric value specifying the minimum score difference
#'   required for a cell to be considered "Positive". This ensures biological
#'   relevance by filtering out weak associations. Default: 0.7.
#' @param verbose Logical, whether to print messages.
#'
#' @return
#' The input `pred_dt` with three additional columns:
#' \itemize{
#'   \item `diff` - Numeric vector of score differences between the two conditions
#'   \item `label` - Character vector with cell classifications: "Positive" or "Other"
#' }
#' The function also provides detailed console output about the classification
#' process and results.
#'
#' @details
#' This function implements a sophisticated approach for binary cell classification
#' that adapts to the underlying distribution of prediction score differences while
#' enforcing a minimum threshold for biological significance:
#'
#' ## Classification Strategies with Minimum Threshold:
#' - **Non-normal distributions (p-value < 0.05)**: Uses quantile-based selection
#'   where the top `select_fraction` of cells by score difference are classified
#'   as "Positive", but only if they exceed `min_threshold`
#' - **Normal distributions (p-value >= 0.05)**: Uses normal distribution
#'   quantiles to determine the classification threshold, adjusted upward if
#'   necessary to meet the minimum threshold requirement
#'
#' ## Supported Normality Tests:
#' - **Jarque-Bera**: Tests for skewness and kurtosis deviations from normality
#' - **D'Agostino**: Extended normality test focusing on skewness
#' - **Kolmogorov-Smirnov**: Non-parametric test comparing empirical distribution
#'   to normal distribution
#'
#' @note
#' The minimum threshold parameter (`min_threshold`) helps prevent over-
#' interpretation of weak phenotypic associations and ensures that classified
#' cells show substantial differences between conditions. The function provides
#' comprehensive feedback about threshold adjustments and final classification
#' statistics.
#'
#' @examples
#' \dontrun{
#' # Create example prediction data
#' pred_data <- data.table(
#'   condition_A = runif(1000),
#'   condition_B = runif(1000)
#' )
#'
#' # Classify cells using D'Agostino test with minimum threshold
#' result <- LabelBinaryCells(
#'   pred_dt = pred_data,
#'   pheno_colnames = c("condition_A", "condition_B"),
#'   select_fraction = 0.1,
#'   test_method = "d'agostino",
#'   min_threshold = 0.7
#' )
#'
#' # Check classification results
#' table(result$label)
#'
#' # View the actual fraction of positive cells
#' prop.table(table(result$label))
#' }
#'
#' @seealso
#' [jb.test.modified()], [dagostino.test()], [ks.test()] for the underlying
#' normality tests used in the classification process.
#'
#' @export
#' @family DEGAS
#'
LabelBinaryCells <- function(
    pred_dt,
    pheno_colnames,
    select_fraction,
    test_method,
    min_threshold = 0.7, # Added minimum threshold parameter
    verbose = SigBridgeRUtils::getFuncOption("verbose")
) {
    chk::chk_length(pheno_colnames, 2)
    # Try to find the reference group in the column names
    # If not found, use the second column as the reference group
    ctrl_col <- grep(
        "[Nn]on|[Nn]ormal|[Cc]ontrol|[Rr]ef|[Cc]trl|0",
        pheno_colnames,
        value = TRUE
    )
    if (nchar(ctrl_col) == 0) {
        if (verbose) {
            cli::cli_alert_info(
                "Using {.val {pheno_colnames[2]}} as reference group"
            )
        }

        pred_dt[, "diff" := .SD[[pheno_colnames[1]]] - .SD[[pheno_colnames[2]]]]
    } else {
        if (verbose) {
            cli::cli_alert_info(
                "Using {.val {ctrl_col}} as reference group"
            )
        }

        pred_dt[,
            "diff" := .SD[[setdiff(pheno_colnames, ctrl_col)]] -
                .SD[[ctrl_col]]
        ]
    }

    # Calculate difference using data.table operations

    # Perform normality test with purrr pattern matching
    normality_test_pval <- switch(
        test_method,
        "jarque-bera" = jb.test.modified(pred_dt$diff)$p.value,
        "d'agostino" = dagostino.test(pred_dt$diff)$p.value[3],
        "kolmogorov-smirnov" = stats::ks.test(pred_dt$diff, "pnorm")$p.value
    )

    # Apply labeling based on normality test result
    pred_dt[,
        "label" := {
            if (normality_test_pval < 0.05) {
                # Use quantile-based selection for non-normal distributions
                n_positive <- ceiling(.N * select_fraction)
                sorted_indices <- order(diff, decreasing = TRUE)
                positive_positions <- sorted_indices[seq_len(n_positive)]

                # Calculate original threshold
                original_thresh <- round(
                    diff[positive_positions[n_positive]],
                    4
                )

                # Apply minimum threshold constraint
                if (original_thresh < min_threshold) {
                    # Filter positions to only include cells above min_threshold
                    above_min_thresh <- which(diff > min_threshold)

                    # Take intersection with top fraction positions
                    valid_positions <- intersect(
                        positive_positions,
                        above_min_thresh
                    )

                    # Update actual number of positive cells
                    actual_n_positive <- length(valid_positions)
                    actual_thresh <- min_threshold

                    if (verbose) {
                        cli::cli_alert_info(
                            "Original threshold {.val {original_thresh}} below minimum {.val {min_threshold}}, using {.val {actual_thresh}} instead"
                        )
                    }
                } else {
                    valid_positions <- positive_positions
                    actual_thresh <- original_thresh
                    actual_n_positive <- n_positive
                }

                if (verbose) {
                    cli::cli_alert_info(
                        "Scores over {.val {actual_thresh}} are considered `Positive`."
                    )
                }

                # Create labels using vectorized operations
                labels <- rep("Other", .N)
                labels[valid_positions] <- "Positive"
                labels
            } else {
                # Use normal distribution-based selection
                mean_val <- colMeans(as.matrix(diff))
                sd_val <- SigBridgeRUtils::colSds(as.matrix(diff))
                quantile_val <- stats::qnorm(
                    select_fraction,
                    mean = mean_val,
                    sd = sd_val,
                    lower.tail = FALSE
                )

                # Apply minimum threshold constraint
                actual_thresh <- max(quantile_val, min_threshold)

                if (quantile_val < min_threshold && verbose) {
                    cli::cli_alert_info(
                        "Original threshold {.val {round(quantile_val, 4)}} below minimum {.val {min_threshold}}, using {.val {actual_thresh}} instead"
                    )
                }
                if (verbose) {
                    cli::cli_alert_info(
                        "Scores over {.val {round(actual_thresh, 4)}} are considered `Positive`."
                    )
                }

                # Use data.table's fast ifelse
                data.table::fifelse(diff > actual_thresh, "Positive", "Other")
            }
        }
    ]

    # Additional validation: count actual positive cells
    positive_count <- pred_dt[label == "Positive", .N]
    total_count <- nrow(pred_dt)
    actual_fraction <- round(positive_count / total_count, 4)
    if (verbose) {
        cli::cli_alert_success(
            "Labeled {.val {positive_count}} cells as Positive ({.val {actual_fraction * 100}}% of total)."
        )
    }

    pred_dt
}

#' @title Label Continuous Phenotype Cells Using MAD Testing
#'
#' @description
#' Identifies phenotype-associated cells in continuous prediction data using
#' Median Absolute Deviation (MAD) testing across multiple phenotypic dimensions.
#'
#' @param pred_dt A data.table containing prediction scores for multiple
#'   phenotypic conditions. Must contain a 'cell_id' column and one or more
#'   columns with prediction scores for different phenotypes.
#' @param verbose Logical, whether to print messages.
#'
#' @return
#' The input `pred_dt` with an additional column:
#' \itemize{
#'   \item `label` - Character vector with cell classifications: "Positive"
#'     (significant in at least one phenotype) or "Other"
#' }
#'
#'
#' @note
#' The function assumes that all columns except 'cell_id' contain prediction
#' scores for different phenotypes. It provides progress information and
#' warnings for edge cases like empty results.
#'
#' @examples
#' \dontrun{
#' # Create example prediction data with multiple phenotypes
#' pred_data <- data.table(
#'   cell_id = paste0("cell_", 1:1000),
#'   phenotype_A = rnorm(1000),
#'   phenotype_B = rexp(1000),
#'   phenotype_C = runif(1000)
#' )
#'
#' # Identify phenotype-associated cells
#' result <- LabelContinuousCells(pred_data)
#'
#' # Check classification results
#' table(result$label)
#'
#' # View the proportion of positive cells
#' prop.table(table(result$label))
#' }
#'
#' @seealso
#' [mad.test()] for the underlying statistical test used for phenotype
#' significance assessment.
#'
#' @export
#' @family DEGAS
#'
LabelContinuousCells <- function(
    pred_dt,
    verbose = SigBridgeRUtils::getFuncOption("verbose")
) {
    if (verbose) {
        ts_cli$cli_alert_info(
            "Searching for various phenotype-associated cells..."
        )
    }

    # Use matrix operations for efficient MAD testing across predictions
    label_cols <- setdiff(names(pred_dt), "cell_id")

    mad_results <- purrr::map_dfr(
        label_cols,
        function(col) {
            vals <- as.numeric(pred_dt[[col]])
            mad_test <- mad.test(vals)
            data.table::data.table(
                column = col,
                p_value = mad_test$p.value,
                is_positive = mad_test$p.value < 0.05
            )
        }
    )
    if (is.null(mad_results)) {
        cli::cli_warn("Empty MAD test results returned.")
        return(pred_dt)
    }

    # Apply labels based on MAD test results
    positive_cols <- mad_results[`p_value` < 0.05, `column`]
    if (length(positive_cols) > 0) {
        pred_dt[,
            "label" := ifelse(rowSums(.SD) > 0, "Positive", "Other"),
            .SDcols = positive_cols
        ]
    } else {
        pred_dt[, "label" := "Other"]
    }

    pred_dt
}

#' @title Label Survival-Associated Phenotype Cells Based on Hazard Scores
#'
#' @description
#' Classifies cells into survival-associated phenotype groups ("Positive" vs "Other")
#' based on hazard scores using statistical distribution analysis. This function
#' identifies cells with significantly elevated hazard scores that may be
#' associated with survival outcomes, employing adaptive thresholding based on
#' distribution characteristics.
#'
#' @param pred_dt A data.table containing hazard scores for cells. Must contain
#'   a column named 'Hazard' with numeric hazard scores.
#' @param select_fraction Numeric value between 0 and 1 specifying the target
#'   fraction of cells to classify as "Positive". The actual fraction may be
#'   adjusted based on distribution characteristics and minimum threshold
#'   constraints.
#' @param test_method Character string specifying the statistical test to use
#'   for normality assessment of hazard scores. One of: `"jarque-bera"`,
#'   `"d'agostino"`, `"kolmogorov-smirnov"`.
#' @param min_threshold minimal threshold for hazard score
#' @param verbose Logical, whether to print messages.
#'
#' @return
#' The input `pred_dt` with an additional column:
#' \itemize{
#'   \item `label` - Character vector with cell classifications: "Positive"
#'     (high hazard cells) or "Other"
#' }
#'
#' @details
#' ## Classification Strategies:
#' - **Non-normal distributions (p-value < 0.05)**: Uses quantile-based selection
#'   where the top `select_fraction` of cells by hazard score are classified as
#'   "Positive", with minimum threshold constraints
#' - **Normal distributions (p-value ≥ 0.05)**: Uses normal distribution
#'   quantiles to determine the classification threshold, adjusted to meet
#'   minimum requirements
#'
#' ## Supported Normality Tests:
#' - **Jarque-Bera**: Tests for skewness and kurtosis deviations from normality
#' - **D'Agostino**: Extended normality test focusing on skewness
#' - **Kolmogorov-Smirnov**: Non-parametric test comparing empirical distribution
#'   to normal distribution
#'
#' @note
#' The function assumes the input data.table contains a column named 'Hazard'
#' with numeric values representing hazard scores from upstream analysis.
#' The minimum threshold is internally defined to ensure biological relevance
#' of the identified cell populations.
#'
#' @examples
#' \dontrun{
#' # Create example hazard score data
#' hazard_data <- data.table(
#'   cell_id = paste0("cell_", 1:1000),
#'   Hazard = rexp(1000, rate = 2)  # Simulated hazard scores
#' )
#'
#' # Identify survival-associated cells
#' result <- LabelSurvivalCells(
#'   pred_dt = hazard_data,
#'   select_fraction = 0.1,
#'   test_method = "jarque-bera"
#' )
#'
#' # Check classification results
#' table(result$label)
#'
#' # Analyze the hazard scores of positive cells
#' summary(result[label == "Positive", Hazard])
#' }
#'
#' @seealso
#' [jb.test.modified()], [dagostino.test()], [ks.test()] for the underlying
#' normality tests used in the classification process.
#'
#' @export
#' @family DEGAS
#'
LabelSurvivalCells <- function(
    pred_dt,
    select_fraction,
    test_method,
    min_threshold = 0.7, # Added minimum threshold parameter
    verbose = SigBridgeRUtils::getFuncOption("verbose")
) {
    if (verbose) {
        ts_cli$cli_alert_info("Searching for survival-associated cells...")
    }

    pred_vec <- pred_dt[["Hazard"]]
    normality_test_pval <- switch(
        test_method,
        "jarque-bera" = jb.test.modified(pred_vec)$p.value,
        "d'agostino" = dagostino.test(pred_vec)$p.value[3],
        "kolmogorov-smirnov" = stats::ks.test(pred_vec, "pnorm")$p.value
    )

    pred_dt[,
        "label" := {
            if (normality_test_pval < 0.05) {
                # Use quantile-based selection for non-normal distributions
                n_positive <- ceiling(.N * select_fraction)
                sorted_indices <- order(`Hazard`, decreasing = TRUE)
                positive_positions <- sorted_indices[seq_len(n_positive)]

                # Calculate original threshold
                original_thresh <- round(
                    `Hazard`[positive_positions[n_positive]],
                    4
                )

                # Apply minimum threshold constraint
                if (original_thresh < min_threshold) {
                    # Filter positions to only include cells above min_threshold
                    above_min_thresh <- which(`Hazard` > min_threshold)

                    # Take intersection with top fraction positions
                    valid_positions <- intersect(
                        positive_positions,
                        above_min_thresh
                    )

                    # * Update actual number of positive cells
                    # actual_n_positive <- length(valid_positions)
                    actual_thresh <- min_threshold

                    if (verbose) {
                        cli::cli_alert_info(
                            "Original threshold {.val {original_thresh}} below minimum {.val {min_threshold}}, using {.val {actual_thresh}} instead"
                        )
                    }
                } else {
                    valid_positions <- positive_positions
                    actual_thresh <- original_thresh
                    # actual_n_positive <- n_positive
                }

                if (verbose) {
                    cli::cli_alert_info(
                        "Scores over {.val {actual_thresh}} are considered `Positive`."
                    )
                }

                # Create labels using vectorized operations
                labels <- rep("Other", .N)
                labels[positive_positions] <- "Positive"
                labels
            } else {
                # Use normal distribution-based selection
                mean_val <- colMeans(`Hazard`)
                sd_val <- SigBridgeRUtils::colSds(`Hazard`)
                quantile_val <- stats::qnorm(
                    select_fraction,
                    mean = mean_val,
                    sd = sd_val,
                    lower.tail = FALSE
                )
                # Apply minimum threshold constraint
                actual_thresh <- max(quantile_val, min_threshold)

                if (quantile_val < min_threshold && verbose) {
                    cli::cli_alert_info(
                        "Original threshold {.val {round(quantile_val, 4)}} below minimum {.val {min_threshold}}, using {.val {actual_thresh}} instead"
                    )
                }

                if (verbose) {
                    cli::cli_alert_info(
                        "Scores over {.val {round(actual_thresh, 4)}} are considered `Positive`."
                    )
                }

                data.table::fifelse(
                    `Hazard` > actual_thresh,
                    "Positive",
                    "Other"
                )
            }
        }
    ]
}
