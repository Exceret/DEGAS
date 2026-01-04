#' @keywords internal
read_py_as_chr <- function(file_name) {
  path <- system.file(
    "DEGAS_tools",
    file_name,
    package = "DEGAS",
    mustWork = TRUE
  )
  paste(readLines(path, warn = FALSE), collapse = "\n")
}
