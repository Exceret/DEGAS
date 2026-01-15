#' @title Optimized Bootstrap Aggregation for Cross-Condition Multi-Task Learning
#'
#' @export
#' @family DEGAS
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
runCCMTLBag.optimized <- function(
  verbose = SigBridgeRUtils::getFuncOption("verbose") %||% TRUE,
  ...
) {
  dots <- rlang::list2(...)
  Bagdepth <- dots$DEGAS.bag_depth %||% 5L
  DEGAS.seed <- dots$DEGAS.seed %||% 123L
  FFdepth <- dots$DEGAS.ff_depth %||% 5L
  architecture <- dots$DEGAS.architecture %||% "DenseNet"
  model_type <- dots$DEGAS.model_type

  if (verbose) {
    ts_cli$cli_alert_info(
      "{FFdepth}-layer {architecture} {model_type} DEGAS model"
    )
  }


  purrr::map(
    seq_len(Bagdepth),
    function(i) {
      DEGAS.seed_i <- DEGAS.seed + (i - 1)

      dots$DEGAS.seed <- DEGAS.seed_i

      if (verbose) {
        ts_cli$cli_alert_info("Training progress: {i}/{Bagdepth}...")
      }

      result <- rlang::exec(runCCMTL.optimized, verbose = verbose, !!!dots, )
      class(result) <- "ccModel"

      result
    },
    .progress = verbose
  )
}
