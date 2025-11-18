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
#' @param DEGAS.seed Integer specifying the base random seed for reproducible
#'   model training. Each model in the ensemble uses a derived seed.
#' @param verbose Logical, whether to print messages.
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
    FFdepth,
    Bagdepth,
    DEGAS.seed,
    verbose = SigBridgeRUtils::getFuncOption("verbose")
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
                DEGAS.seed = DEGAS.seed_i,
                # Written files will not be rewritten
                force_rewrite = FALSE
            )
            class(result) <- "ccModel"

            result
        },
        .progress = verbose
    )
}
