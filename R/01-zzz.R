# ? Package startup messages
.onLoad <- function(libname, pkgname) {
  # Declare Python dependencies for reticulate
  reticulate::py_require(c("tensorflow", "numpy"))
  invisible()
}

.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion(pkgname)

  msg <- cli::cli_fmt(cli::cli_alert_success(
    "{.pkg {pkgname}} v{pkg_version} loaded"
  ))
  packageStartupMessage(msg)
  invisible()
}

ts_cli <- SigBridgeRUtils::CreateTimeStampCliEnv()

py <- reticulate::py

#'@importFrom data.table `%chin%` %chin%
NULL
