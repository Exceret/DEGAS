#' @title Optimized Input File Writing for DEGAS Models
#'
#' @description
#' Efficiently writes input data files for DEGAS model training using optimized
#' data handling and fast I/O operations. This function converts various data
#' types to efficient CSV format using data.table's fwrite for rapid file
#' operations with comprehensive error handling.
#'
#' @param scExp A matrix, data frame, or Matrix object containing single-cell
#'   expression data. Rows typically represent genes and columns represent cells.
#' @param scLab A matrix, or data frame containing single-cell labels
#'   corresponding to the expression data. Can be NULL if no labels are available.
#' @param patExp A data frame, or Matrix object containing patient-level
#'   expression data. Rows typically represent genes and columns represent patients.
#' @param patLab A matrix, or data frame containing patient-level labels
#'   corresponding to the patient expression data. Can be NULL if no labels are available.
#' @param tmpDir Character string specifying the directory path where input
#'   files will be written. The directory will be created if it doesn't exist.
#'
#' @return
#' Invisibly returns TRUE if all files are successfully written. If any error
#' occurs during file writing, the function will abort with an informative error message.
#'
#' @details
#' This function provides an optimized pipeline for writing input files required
#' by DEGAS models with the following features:
#'
#' ## File Output:
#' The function creates four CSV files in the specified temporary directory:
#' - `scExp.csv`: Single-cell expression data
#' - `scLab.csv`: Single-cell labels (if provided)
#' - `patExp.csv`: Patient-level expression data
#'
#'
#' @note
#' The function uses comma-separated values (CSV) format without row names to
#' ensure compatibility with Python-based DEGAS training scripts. All input
#' data is converted to dense format during writing, so ensure sufficient
#' memory is available for large datasets.
#'
#' @examples
#' \dontrun{
#' # Write input files for DEGAS training
#' writeInputFiles.optimized(
#'   scExp = single_cell_expression,
#'   scLab = single_cell_labels,
#'   patExp = patient_expression,
#'   patLab = patient_labels,
#'   tmpDir = "/tmp/degas_input"
#' )
#' }
#'
#' @seealso
#' [data.table::fwrite()] for the underlying fast writing implementation,
#' [purrr::safely()] for the error handling mechanism.
#'
#' @export
#' @family DEGAS
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
writeInputFiles.optimized <- function(
    scExp,
    scLab = NULL,
    patExp,
    patLab = NULL,
    tmpDir
) {
    if (!dir.exists(tmpDir)) {
        dir.create(tmpDir, recursive = TRUE)
    }

    file_writing_ops <- purrr::safely(
        ~ {
            # Convert scExp to data.table for fast writing
            if (inherits(scExp, "matrix") || inherits(scExp, "Matrix")) {
                scExp_dt <- data.table::as.data.table(as.matrix(scExp))
            } else {
                scExp_dt <- data.table::as.data.table(scExp)
            }

            # Write scExp using data.table's fwrite (much faster than write.table)
            data.table::fwrite(
                scExp_dt,
                file = file.path(tmpDir, 'scExp.csv'),
                row.names = FALSE,
                sep = ',',
                showProgress = FALSE
            )

            # Process scLab if not NULL
            if (!is.null(scLab)) {
                if (inherits(scLab, "matrix") || inherits(scLab, "Matrix")) {
                    scLab_dt <- data.table::as.data.table(as.matrix(scLab))
                } else {
                    scLab_dt <- data.table::as.data.table(scLab)
                }
                scLab_dt[, names(scLab_dt) := lapply(.SD, as.integer)] # Convert to integer from boolean

                data.table::fwrite(
                    scLab_dt,
                    file = file.path(tmpDir, 'scLab.csv'),
                    row.names = FALSE,
                    sep = ',',
                    showProgress = FALSE
                )
            }

            # Convert patExp to data.table
            if (inherits(patExp, "matrix") || inherits(patExp, "Matrix")) {
                patExp_dt <- data.table::as.data.table(as.matrix(patExp))
            } else {
                patExp_dt <- data.table::as.data.table(patExp)
            }

            # Write patExp
            data.table::fwrite(
                patExp_dt,
                file = file.path(tmpDir, 'patExp.csv'),
                row.names = FALSE,
                sep = ',',
                showProgress = FALSE
            )

            # Process patLab if not NULL
            if (!is.null(patLab)) {
                if (inherits(patLab, "matrix") || inherits(patLab, "Matrix")) {
                    patLab_dt <- data.table::as.data.table(as.matrix(patLab))
                } else {
                    patLab_dt <- data.table::as.data.table(patLab)
                }
                patLab_dt[, names(patLab_dt) := lapply(.SD, as.integer)] # Convert to integer from boolean

                data.table::fwrite(
                    patLab_dt,
                    file = file.path(tmpDir, 'patLab.csv'),
                    row.names = FALSE,
                    sep = ',',
                    showProgress = FALSE
                )
            }
        }
    )

    # Execute file writing operations with error handling
    result <- file_writing_ops()

    if (!is.null(result$error)) {
        cli::cli_abort("Failed to write input files: ", result$error$message)
    }

    invisible(TRUE)
}
