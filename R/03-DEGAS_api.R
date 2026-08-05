#' Train DEGAS model
#'
#' @param sc_matrix Single-cell expression matrix, cells x genes.
#' @param bulk_matrix Bulk expression matrix, samples x genes.
#' @param bulk_labels Cox labels or classification labels.
#'   Cox: data.frame with columns `time` and `status`.
#'   Classification: named list/vector or dict-like object, names are bulk sample names.
#' @param sc_gene_names Gene names for columns of sc_matrix.
#' @param bulk_gene_names Gene names for columns of bulk_matrix.
#' @param bulk_sample_names Sample names for rows of bulk_matrix.
#' @param cell_names Cell names for rows of sc_matrix.
#' @param cell_labels Optional cell type labels.
#' @param patient_task "auto", "cox", or "classification".
#' @param convert Whether reticulate converts Python objects to R.
#'   For model training this should usually be FALSE so the Python model object is preserved.
#' @param ... Additional arguments passed to Python DEGASTensorFlow.
#'
#' @return Python DEGASTensorFlow model object.
#' @export
train_degas <- function(
  sc_matrix,
  bulk_matrix,
  bulk_labels,
  sc_gene_names,
  bulk_gene_names,
  bulk_sample_names = NULL,
  cell_names = NULL,
  cell_labels = NULL,
  patient_task = c("auto", "cox", "classification"),
  convert = FALSE,
  ...
) {
  patient_task <- match.arg(patient_task)

  py <- degas_py(convert = convert)

  py$train_degas(
    sc_matrix = sc_matrix,
    bulk_matrix = bulk_matrix,
    bulk_labels = bulk_labels,
    sc_gene_names = sc_gene_names,
    bulk_gene_names = bulk_gene_names,
    bulk_sample_names = bulk_sample_names,
    cell_names = cell_names,
    cell_labels = cell_labels,
    patient_task = patient_task,
    ...
  )
}


#' Predict cell-level DEGAS scores
#'
#' @param model A Python DEGASTensorFlow model returned by train_degas.
#' @param sc_matrix Single-cell expression matrix, cells x genes.
#' @param sc_gene_names Gene names for columns of sc_matrix.
#' @param cell_names Optional cell names.
#' @param batch_size Prediction batch size.
#' @param convert Whether to convert Python pandas.DataFrame to R data.frame.
#'
#' @return data.frame if convert=TRUE, otherwise Python pandas.DataFrame.
#' @export
predict_degas <- function(
  model,
  sc_matrix,
  sc_gene_names,
  cell_names = NULL,
  batch_size = 2048L,
  convert = TRUE
) {
  py <- degas_py(convert = convert)

  py$predict_degas(
    model = model,
    sc_matrix = sc_matrix,
    sc_gene_names = sc_gene_names,
    cell_names = cell_names,
    batch_size = as.integer(batch_size)
  )
}


#' Run DEGAS training and prediction in one call
#'
#' @param sc_matrix Single-cell expression matrix, cells x genes.
#' @param bulk_matrix Bulk expression matrix, samples x genes.
#' @param bulk_labels Cox labels or classification labels.
#' @param sc_gene_names Gene names for columns of sc_matrix.
#' @param bulk_gene_names Gene names for columns of bulk_matrix.
#' @param bulk_sample_names Sample names for rows of bulk_matrix.
#' @param cell_names Cell names for rows of sc_matrix.
#' @param cell_labels Optional cell type labels.
#' @param patient_task "auto", "cox", or "classification".
#' @param convert Whether to convert output to R data.frame.
#' @param ... Additional parameters passed to Python backend.
#'
#' @return Cell-level DEGAS prediction scores.
#' @export
run_degas <- function(
  sc_matrix,
  bulk_matrix,
  bulk_labels,
  sc_gene_names,
  bulk_gene_names,
  bulk_sample_names = NULL,
  cell_names = NULL,
  cell_labels = NULL,
  patient_task = c("auto", "cox", "classification"),
  convert = TRUE,
  ...
) {
  patient_task <- match.arg(patient_task)

  py <- degas_py(convert = convert)

  py$run_degas(
    sc_matrix = sc_matrix,
    bulk_matrix = bulk_matrix,
    bulk_labels = bulk_labels,
    sc_gene_names = sc_gene_names,
    bulk_gene_names = bulk_gene_names,
    bulk_sample_names = bulk_sample_names,
    cell_names = cell_names,
    cell_labels = cell_labels,
    patient_task = patient_task,
    ...
  )
}
