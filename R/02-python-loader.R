# ? Lazy loader for the bundled Python DEGAS backend
#
# The Python module under inst/DEGAS is imported lazily on first use so that
# loading or documenting the R package never requires a working Python
# environment (devtools::document()/load_all() runs .onLoad(), which must not
# touch Python). Use degas_py() inside package functions instead of importing
# in .onLoad().

.degas_py_cache <- new.env(parent = emptyenv())

#' Load bundled Python DEGAS module
#'
#' This function lazily imports the Python module under `inst/DEGAS`.
#' Users normally do not need to call this directly.
#'
#' @param convert Whether reticulate should convert Python objects to R objects.
#' @param reload Force re-importing the Python module.
#'
#' @return Python module object.
#' @export
degas_py <- function(convert = TRUE, reload = FALSE) {
  if (reload) {
    degas_py_reset()
  }
  cached <- exists("module", envir = .degas_py_cache, inherits = FALSE) &&
    identical(
      get("convert", envir = .degas_py_cache, inherits = FALSE),
      convert
    )
  if (cached) {
    return(get("module", envir = .degas_py_cache))
  }

  py <- reticulate::import_from_path(
    module = "DEGAS",
    path = system.file(package = "DEGAS")
  )
  assign("module", py, envir = .degas_py_cache)
  assign("convert", convert, envir = .degas_py_cache)
  py
}

#' Check whether bundled Python DEGAS backend can be loaded
#'
#' @return TRUE/FALSE.
#' @export
degas_py_available <- function() {
  if (!reticulate::py_available(initialize = FALSE)) {
    return(FALSE)
  }
  tryCatch(
    {
      degas_py()
      TRUE
    },
    error = function(e) FALSE
  )
}

#' Reset cached Python DEGAS module
#'
#' This only clears the R-side cached module reference.
#'
#' @return Invisible NULL.
#' @export
degas_py_reset <- function() {
  if (exists("module", envir = .degas_py_cache, inherits = FALSE)) {
    rm("module", envir = .degas_py_cache)
    rm("convert", envir = .degas_py_cache)
  }
  invisible(NULL)
}
