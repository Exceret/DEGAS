#' @title Optimized Bootstrap Aggregation for Cross-Condition Multi-Task Learning
#'
#' @description
#' Performs bootstrap aggregated training of multiple CCMTL models to enhance
#' robustness and reduce variance in predictions. This function trains an
#' ensemble of models with different random seeds and aggregates the results.
#'
#' @param scExp A matrix or data frame containing single-cell expression data
#'   for model training.
#' @param scLab A matrix containing single-cell labels corresponding
#'   to the expression data.
#' @param patExp A matrix or data frame containing patient-level expression data
#'   for multi-task learning.
#' @param patLab A matrix containing patient-level labels corresponding
#'   to the patient expression data.
#' @param tmpDir Character string specifying the temporary directory path for
#'   storing intermediate files and model outputs.
#' @param model_type Character string specifying the type of model to train.
#'   Should match available DEGAS model types.
#' @param architecture Character string specifying the neural network architecture.
#'   One of: "DenseNet", "Standard".
#' @param FFdepth Integer specifying the number of layers in the feed-forward
#'   network architecture.
#' @param Bagdepth Integer specifying the number of bootstrap models to train
#'   in the ensemble.
#' @param DEGAS.pyloc Character string specifying the path to the Python
#'   executable to use for model training.
#' @param DEGAS.toolsPath Character string specifying the path to the DEGAS
#'   Python tools directory.
#' @param DEGAS.train_steps Integer specifying the number of training steps.
#'   Default: 2000.
#' @param DEGAS.scbatch_sz Integer specifying the single-cell batch size.
#'   Default: 200.
#' @param DEGAS.patbatch_sz Integer specifying the patient batch size.
#'   Default: 50.
#' @param DEGAS.hidden_feats Integer specifying the number of hidden features.
#'   Default: 50.
#' @param DEGAS.do_prc Numeric specifying the dropout keep probability (0-1).
#'   Default: 0.5.
#' @param DEGAS.lambda1 Numeric specifying the L2 regularization term.
#'   Default: 3.0.
#' @param DEGAS.lambda2 Numeric specifying the patient loss term.
#'   Default: 3.0.
#' @param DEGAS.lambda3 Numeric specifying the MMD loss term.
#'   Default: 3.0.
#' @param DEGAS.seed Integer specifying the base random seed for reproducible
#'   model training. Each model in the ensemble uses a derived seed.
#' @param verbose Logical, whether to print messages.
#' @param force_rewrite rewrite input files
#' @param ... unused
#'
#' @return
#' Returns a list of trained CCMTL model objects from the bootstrap aggregation
#' process. The list contains successful model results with proper error handling
#' for failed training attempts.
#'
#' @details
#' This function implements bootstrap aggregated training (bagging) for CCMTL
#' models with the following features:
#'
#' ## Ensemble Training:
#' - Trains multiple models with different random seeds derived from the base seed
#' - Uses parallel-safe file management to avoid I/O conflicts
#' - Implements comprehensive error handling to continue training even if
#'   individual models fail
#'
#' ## Error Handling:
#' - Continues training even if individual models fail
#' - Returns only successfully trained models
#' - Provides progress feedback for long-running ensemble training
#'
#' @note
#' The bootstrap aggregation process can be computationally intensive, especially
#' for large datasets or deep architectures. The function creates derived seeds
#' for each model (base seed + model index) to ensure reproducibility while
#' maintaining diversity in the ensemble.
#'
#' @examples
#' \dontrun{
#' # Train an ensemble of 10 CCMTL models
#' ensemble_models <- runCCMTLBag.optimized(
#'   scExp = sc_expression,
#'   scLab = sc_labels,
#'   patExp = patient_expression,
#'   patLab = patient_labels,
#'   tmpDir = "/tmp/degas_models",
#'   model_type = "classification",
#'   architecture = "DenseNet",
#'   FFdepth = 3,
#'   Bagdepth = 10,
#'   DEGAS.pyloc = "python3",
#'   DEGAS.toolsPath = "/path/to/tools/",
#'   DEGAS.seed = 42
#' )
#'
#' # Access individual models from the ensemble
#' first_model <- ensemble_models[[1]]
#' }
#'
#' @seealso
#' [runCCMTL.optimized] for single model training,
#' [purrr::map()] for the iterative execution pattern.
#'
#' @export
#' @family DEGAS
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
runCCMTLBag.optimized <- function(
  scExp,
  scLab,
  patExp,
  patLab,
  tmpDir,
  model_type,
  architecture,
  FFdepth = 3L,
  Bagdepth = 5L,
  DEGAS.pyloc,
  DEGAS.toolsPath,
  DEGAS.train_steps = 2000L,
  DEGAS.scbatch_sz = 200L,
  DEGAS.patbatch_sz = 50L,
  DEGAS.hidden_feats = 50L,
  DEGAS.do_prc = 0.5,
  DEGAS.lambda1 = 3.0,
  DEGAS.lambda2 = 3.0,
  DEGAS.lambda3 = 3.0,
  DEGAS.seed,
  verbose = SigBridgeRUtils::getFuncOption("verbose") %||% TRUE,
  force_rewrite = FALSE,
  ...
) {
  if (verbose) {
    ts_cli$cli_alert_info(
      "{FFdepth}-layer {architecture} {model_type} DEGAS model"
    )
  }

  if (!dir.exists(tmpDir)) {
    dir.create(tmpDir, recursive = TRUE)
  }
  # Write files once at the beginning
  writeInputFiles.optimized(
    scExp = scExp,
    scLab = scLab,
    patExp = patExp,
    patLab = patLab,
    tmpDir = tmpDir
  )

  # Check Python with error handling
  py_check <- processx::run(command = DEGAS.pyloc, args = "--version")
  if (!is.null(py_check$error)) {
    cli::cli_abort("Python check failed: ", py_check$error$message)
  } else if (verbose) {
    ts_cli$cli_alert_info(
      "Python check passed, using {py_check$stdout}"
    )
  }

  if (mirai::daemons_set()) {
    ccmtl_fn <- purrr::in_parallel(
      function(i) {
        DEGAS.seed_i <- DEGAS.seed + (i - 1)

        if (verbose) {
          ts_cli$cli_alert_info("Training progress: {i}/{Bagdepth}...")
        }

        result <- runCCMTL.optimized(
          scExp = scExp,
          scLab = scLab,
          patExp = patExp,
          patLab = patLab,
          tmpDir = tmpDir,
          model_type = model_type,
          architecture = architecture,
          FFdepth = FFdepth,
          DEGAS.pyloc = DEGAS.pyloc,
          DEGAS.toolsPath = DEGAS.toolsPath,
          DEGAS.train_steps = DEGAS.train_steps,
          DEGAS.scbatch_sz = DEGAS.scbatch_sz,
          DEGAS.patbatch_sz = DEGAS.patbatch_sz,
          DEGAS.hidden_feats = DEGAS.hidden_feats,
          DEGAS.do_prc = DEGAS.do_prc,
          DEGAS.lambda1 = DEGAS.lambda1,
          DEGAS.lambda2 = DEGAS.lambda2,
          DEGAS.lambda3 = DEGAS.lambda3,
          DEGAS.seed = DEGAS.seed_i,
          # Written files will not be rewritten
          force_rewrite = force_rewrite
        )
        class(result) <- "ccModel"

        result
      },
      Bagdepth = Bagdepth,
      ts_cli = ts_cli,
      verbose = verbose,
      runCCMTL.optimized = runCCMTL.optimized,
      scExp = scExp,
      scLab = scLab,
      patExp = patExp,
      patLab = patLab,
      tmpDir = tmpDir,
      model_type = model_type,
      architecture = architecture,
      FFdepth = FFdepth,
      DEGAS.pyloc = DEGAS.pyloc,
      DEGAS.toolsPath = DEGAS.toolsPath,
      DEGAS.train_steps = DEGAS.train_steps,
      DEGAS.scbatch_sz = DEGAS.scbatch_sz,
      DEGAS.patbatch_sz = DEGAS.patbatch_sz,
      DEGAS.hidden_feats = DEGAS.hidden_feats,
      DEGAS.do_prc = DEGAS.do_prc,
      DEGAS.lambda1 = DEGAS.lambda1,
      DEGAS.lambda2 = DEGAS.lambda2,
      DEGAS.lambda3 = DEGAS.lambda3,
      DEGAS.seed = DEGAS.seed,
      force_rewrite = force_rewrite,
      makeExec = makeExec,
      makeExec2 = makeExec2
    )

    purrr::map(
      seq_len(Bagdepth),
      ccmtl_fn,
      .progress = verbose
    )
  } else {
    purrr::map(
      seq_len(Bagdepth),
      function(i) {
        DEGAS.seed_i <- DEGAS.seed + (i - 1)

        if (verbose) {
          ts_cli$cli_alert_info("Training progress: {i}/{Bagdepth}...")
        }

        result <- runCCMTL.optimized(
          scExp = scExp,
          scLab = scLab,
          patExp = patExp,
          patLab = patLab,
          tmpDir = tmpDir,
          model_type = model_type,
          architecture = architecture,
          FFdepth = FFdepth,
          DEGAS.pyloc = DEGAS.pyloc,
          DEGAS.toolsPath = DEGAS.toolsPath,
          DEGAS.train_steps = DEGAS.train_steps,
          DEGAS.scbatch_sz = DEGAS.scbatch_sz,
          DEGAS.patbatch_sz = DEGAS.patbatch_sz,
          DEGAS.hidden_feats = DEGAS.hidden_feats,
          DEGAS.do_prc = DEGAS.do_prc,
          DEGAS.lambda1 = DEGAS.lambda1,
          DEGAS.lambda2 = DEGAS.lambda2,
          DEGAS.lambda3 = DEGAS.lambda3,
          DEGAS.seed = DEGAS.seed_i,
          # Written files will not be rewritten
          force_rewrite = force_rewrite
        )
        class(result) <- "ccModel"

        result
      },
      .progress = verbose
    )
  }
}
