#' @title Optimized Bootstrap Aggregation for Cross-Condition Multi-Task Learning
#'
#' @export
#' @family DEGAS
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
runCCMTLBag.optimized <- function(
  verbose = SigBridgeRUtils::getFuncOption("verbose") %||% TRUE,
  #   scExp = matrix(), # sc data matrix
  #   scLab = matrix(),
  #   patExp = matrix(),
  #   patLab = matrix(),
  #   #   tmpDir,
  #   DEGAS.model_type = c(
  #     'ClassClass',
  #     'ClassCox',
  #     'ClassBlank',
  #     'BlankClass',
  #     'BlankCox'
  #   ),
  #   DEGAS.architecture = c("DenseNet", "Standard"),
  #   DEGAS.ff_depth = 3L,
  #   DEGAS.bag_depth = 5L,
  #   DEGAS.pyloc = ListPyEnv()$python[1],
  #   DEGAS.toolsPath = file.path(.libPaths()[1], "DEGAS/DEGAS_tools/"),
  #   DEGAS.train_steps = 2000L,
  #   DEGAS.scbatch_sz = 200L,
  #   DEGAS.patbatch_sz = 50L,
  #   DEGAS.hidden_feats = 50L,
  #   DEGAS.do_prc = 0.5,
  #   DEGAS.lambda1 = 3.0,
  #   DEGAS.lambda2 = 3.0,
  #   DEGAS.lambda3 = 3.0,
  #   DEGAS.seed = 123L,
  ... # path.data, path.result, etc
) {
  if (verbose) {
    ts_cli$cli_alert_info(
      "{FFdepth}-layer {architecture} {model_type} DEGAS model"
    )
  }

  #   if (!dir.exists(tmpDir)) {
  #     dir.create(tmpDir, recursive = TRUE)
  #   }
  #   # Write files once at the beginning
  #   writeInputFiles.optimized(
  #     scExp = scExp,
  #     scLab = scLab,
  #     patExp = patExp,
  #     patLab = patLab,
  #     tmpDir = tmpDir
  #   )

  # Check Python with error handling
  #   py_check <- processx::run(command = DEGAS.pyloc, args = "--version")
  #   if (!is.null(py_check$error)) {
  #     cli::cli_abort("Python check failed: ", py_check$error$message)
  #   } else if (verbose) {
  #     ts_cli$cli_alert_info(
  #       "Python check passed, using {py_check$stdout}"
  #     )
  #   }

  purrr::map(
    seq_len(Bagdepth),
    function(i) {
      DEGAS.seed_i <- DEGAS.seed + (i - 1)

      if (verbose) {
        ts_cli$cli_alert_info("Training progress: {i}/{Bagdepth}...")
      }

      result <- runCCMTL.optimized(
        verbose = verbose,
        # scExp = scExp, # sc data matrix
        # scLab = scLab,
        # patExp = patExp,
        # patLab = patLab,
        # # tmpDir,
        # DEGAS.model_type = DEGAS.model_type,
        # DEGAS.architecture = DEGAS.architecture,
        # DEGAS.ff_depth = DEGAS.ff_depth,
        # DEGAS.bag_depth = DEGAS.bag_depth,
        # DEGAS.pyloc = DEGAS.pyloc,
        # DEGAS.toolsPath = DEGAS.toolsPath,
        # DEGAS.train_steps = DEGAS.train_steps,
        # DEGAS.scbatch_sz = DEGAS.scbatch_sz,
        # DEGAS.patbatch_sz = DEGAS.patbatch_sz,
        # DEGAS.hidden_feats = DEGAS.hidden_feats,
        # DEGAS.do_prc = DEGAS.do_prc,
        # DEGAS.lambda1 = DEGAS.lambda1,
        # DEGAS.lambda2 = DEGAS.lambda2,
        # DEGAS.lambda3 = DEGAS.lambda3,
        # DEGAS.seed = DEGAS.seed,
        ... # path.data, path.result, etc
      )
      class(result) <- "ccModel"

      result
    },
    .progress = verbose
  )
}
