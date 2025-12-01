#' @title Optimized Cross-Condition Multi-Task Learning Model Training
#'
#' @description
#' An optimized wrapper function for training cross-condition multi-task learning
#' (CCMTL) models in the DEGAS framework. This function handles the complete
#' training pipeline including data preparation, model configuration, and
#' execution with enhanced performance and error handling.
#'
#' @param scExp A matrix or data frame containing single-cell expression data
#'   for model training.
#' @param scLab A  matrix containing single-cell labels corresponding
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
#' @param DEGAS.seed Integer specifying the random seed for reproducible
#'   model training.
#' @param force_rewrite Logical indicating whether to force rewriting of input
#'   files even if they already exist. Default: FALSE.
#' @param verbose Logical, whether to print output messages.
#'
#' @return
#' Returns a trained CCMTL model object that can be used for predictions and
#' further analysis.
#'
#' @details
#' ## Workflow:
#' 1. **File Management**: Efficient handling of temporary directories and
#'    input files with optional forced rewriting
#' 2. **Architecture Configuration**: Supports multiple neural network
#'    architectures (DenseNet, Standard) with customizable depth
#' 3. **Python Environment**: Validates Python availability and executes
#'    training scripts with proper error handling
#' 4. **Model Training**: Executes the DEGAS training process with specified
#'    hyperparameters and architecture choices
#'
#' @note
#' This function requires a properly configured Python environment with DEGAS
#' dependencies installed. The temporary directory (`tmpDir`) should have
#' sufficient disk space for model files and intermediate data.
#'
#' @seealso
#' [runCCMTLBag.optimized()] for bootstrap aggregated model training,
#'
#' @export
#' @family DEGAS
#'
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
runCCMTL.optimized <- function(
    scExp,
    scLab,
    patExp,
    patLab,
    tmpDir,
    model_type,
    architecture,
    FFdepth,
    DEGAS.seed,
    force_rewrite = FALSE,
    verbose = SigBridgeRUtils::getFuncOption("verbose") %||% TRUE
) {
    # Only write files if explicitly requested
    if (force_rewrite) {
        if (dir.exists(tmpDir)) {
            unlink(tmpDir, recursive = TRUE, force = TRUE)
        }
        dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)
        # Write input files
        writeInputFiles.optimized(
            scExp = scExp,
            scLab = scLab,
            patExp = patExp,
            patLab = patLab,
            tmpDir = tmpDir
        )
    }

    # create python files
    if (!architecture %chin% c("DenseNet", "Standard")) {
        cli::cli_abort(c("x" = 'Incorrect architecture argument'))
    } else if (architecture == "DenseNet") {
        makeExec2(
            tmpDir = tmpDir,
            FFdepth = FFdepth,
            model_type = model_type
        )
    } else {
        makeExec(
            tmpDir = tmpDir,
            FFdepth = FFdepth,
            model_type = model_type
        )
    }

    # cmd <- paste0(
    #     DEGAS.pyloc,
    #     " ",
    #     tmpDir,
    #     model_type,
    #     "MTL.py",
    #     paste(
    #         "",
    #         tmpDir,
    #         DEGAS.train_steps,
    #         DEGAS.scbatch_sz,
    #         DEGAS.patbatch_sz,
    #         DEGAS.hidden_feats,
    #         DEGAS.do_prc,
    #         DEGAS.lambda1,
    #         DEGAS.lambda2,
    #         DEGAS.lambda3,
    #         DEGAS.seed
    #     )
    # )
    cmd_args <- as.character(c(
        paste0(tmpDir, model_type, "MTL.py"),
        tmpDir,
        DEGAS.train_steps,
        DEGAS.scbatch_sz,
        DEGAS.patbatch_sz,
        DEGAS.hidden_feats,
        DEGAS.do_prc,
        DEGAS.lambda1,
        DEGAS.lambda2,
        DEGAS.lambda3,
        DEGAS.seed
    ))

    # * Execute system command
    # system(command = cmd) # if processx::run failed, use `system` instead
    result <- processx::run(
        command = DEGAS.pyloc,
        args = cmd_args,
        echo = verbose,
        error_on_status = TRUE
    )

    readOutputFiles.optimized(
        tmpDir = tmpDir,
        model_type = model_type,
        architecture = architecture
    )
}
