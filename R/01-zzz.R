# ? Package startup messages
.onLoad <- function(libname, pkgname) {
  # Declare Python dependencies so reticulate automatically configures a
  # Python environment (installing any missing packages) on first use.
  #
  # NOTE: we deliberately do NOT import the DEGAS Python module here. Doing so
  # triggers a real Python import during devtools::document()/load_all() and
  # fails when the configured Python lacks scipy/tensorflow. The module is
  # imported lazily on first use via degas_py() (see R/02-python-loader.R).
  reticulate::py_require(c(
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "scikit-survival",
    "tensorflow",
    "keras"
  ))

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
