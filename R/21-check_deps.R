# Return operating system

#' @title System and Dependency Checks
#' @name system-checks
#' @description
#' Utility functions for checking operating system and Python dependencies.
NULL

#' @rdname system-checks
#' @description
#' checkOS returns the operating system name.
#'
#' @return Character string of the operating system (e.g., "Windows", "Linux", "Darwin")
#'
#' @examples
#' # Check current operating system
#' os <- checkOS()
#' print(os)
#'
#' @export
checkOS <- function() {
  return(Sys.info()['sysname'])
}

#' @rdname system-checks
#' @description
#' checkForPy checks Python version using the configured Python location.
#'
#' @param pyloc Character string specifying the path to the Python executable.
#'
#' @return System command exit status (0 for success)
#'
#' @examples
#' # Check Python version
#' \dontrun{
#' if (checkForPy("python3") == 0) {
#' print("Python is available")
#' }
#' }
#'
#' @export
checkForPy <- function(pyloc) {
  return(system(paste0(pyloc, " -V")))
}

#' @rdname system-checks
#' @description
#' checkForTF checks if TensorFlow can be imported in Python.
#'
#' @param pyloc Character string specifying the path to the Python executable.
#'
#' @return System command exit status (0 if TensorFlow imports successfully)
#'
#' @examples
#' # Check TensorFlow availability
#' \dontrun{
#' if (checkForTF("python3") == 0) {
#' print("TensorFlow is available")
#' }
#' }
#'
#' @export
checkForTF <- function(pyloc) {
  return(system(paste0(pyloc, " -c 'import tensorflow'")))
}
