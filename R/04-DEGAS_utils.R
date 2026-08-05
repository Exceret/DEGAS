#' Convert DEGAS probabilities to association scores
#'
#' @param probs Numeric vector or matrix.
#' @param convert Whether to convert Python result to R.
#'
#' @return Numeric vector or matrix.
#' @export
toCorrCoeff <- function(probs, convert = TRUE) {
  py <- degas_py(convert = convert)
  py$toCorrCoeff(probs)
}


#' Center DEGAS scores
#'
#' @param x Numeric vector or matrix.
#' @param convert Whether to convert Python result to R.
#'
#' @return Centered vector or matrix.
#' @export
centerFunc <- function(x, convert = TRUE) {
  py <- degas_py(convert = convert)
  py$centerFunc(x)
}


#' kNN smooth DEGAS scores
#'
#' @param probs DEGAS scores, vector or matrix.
#' @param locs Coordinates, for example UMAP or tSNE matrix.
#' @param k Number of nearest neighbors.
#' @param convert Whether to convert Python result to R.
#'
#' @return Smoothed scores.
#' @export
knnSmooth <- function(probs, locs, k = 5L, convert = TRUE) {
  py <- degas_py(convert = convert)

  py$knnSmooth(
    probs = probs,
    locs = locs,
    k = as.integer(k)
  )
}


#' Preprocess count matrix using DEGAS-style transform
#'
#' @param X Count matrix, genes x samples/cells.
#' @param convert Whether to convert Python result to R.
#'
#' @return Processed matrix, samples/cells x genes.
#' @export
preprocessCounts <- function(X, convert = TRUE) {
  py <- degas_py(convert = convert)
  py$preprocessCounts(X)
}
