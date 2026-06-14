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
  verbose = SigBridgeRUtils::getFuncOption("verbose") %||% TRUE,
  scExp = matrix(), # sc data matrix
  scLab = matrix(),
  patExp = matrix(),
  patLab = matrix(),
  #   tmpDir,
  DEGAS.model_type = c(
    "ClassClass",
    "ClassCox",
    "ClassBlank",
    "BlankClass",
    "BlankCox"
  ),
  DEGAS.architecture = c("DenseNet", "Standard"),
  DEGAS.ff_depth = 3L,
  DEGAS.bag_depth = 5L,
  DEGAS.pyloc = ListPyEnv()$python[1],
  DEGAS.toolsPath = file.path(.libPaths()[1], "DEGAS/DEGAS_tools/"),
  DEGAS.train_steps = 2000L,
  DEGAS.scbatch_sz = 200L,
  DEGAS.patbatch_sz = 50L,
  DEGAS.hidden_feats = 50L,
  DEGAS.do_prc = 0.5,
  DEGAS.lambda1 = 3.0,
  DEGAS.lambda2 = 3.0,
  DEGAS.lambda3 = 3.0,
  DEGAS.seed = 123L,
  ... # path.data, path.result, assay etc
) {
  # create python files
  if (!DEGAS.architecture %chin% c("DenseNet", "Standard")) {
    cli::cli_abort(c(
      "x" = "Incorrect architecture argument",
      ">" = "Available architectures: 'DenseNet', 'Standard'"
    ))
  } else {
    full_degas_script <- makeExec(
      #   tmpDir = tmpDir,
      FFdepth = DEGAS.ff_depth,
      model_type = DEGAS.model_type,
      DEGAS.toolsPath = DEGAS.toolsPath,
      architecture = DEGAS.architecture,
      # "Standard", # ! makeExec
      # "DenseNet" # ! makeExec2
    )
  }
  reticulate::py_run_string("import os") # startup
  py <- reticulate::py

  # R matrix -> nparray (set individual globals for Python access)
  py$Xsc <- reticulate::r_to_py(scExp)
  py$Ysc <- reticulate::r_to_py(scLab)
  py$Xpat <- reticulate::r_to_py(patExp)
  py$Ypat <- reticulate::r_to_py(patLab)
  py$train_steps <- reticulate::r_to_py(DEGAS.train_steps)
  py$scbatch_sz <- reticulate::r_to_py(DEGAS.scbatch_sz)
  py$patbatch_sz <- reticulate::r_to_py(DEGAS.patbatch_sz)
  py$hidden_feats <- reticulate::r_to_py(DEGAS.hidden_feats)
  py$do_prc <- reticulate::r_to_py(DEGAS.do_prc)
  py$lambda1 <- reticulate::r_to_py(DEGAS.lambda1)
  py$lambda2 <- reticulate::r_to_py(DEGAS.lambda2)
  py$lambda3 <- reticulate::r_to_py(DEGAS.lambda3)
  py$seed <- reticulate::r_to_py(DEGAS.seed)

  reticulate::py_run_string(full_degas_script)

  # * extract result
  activation <- py$activation # A list

  additional_layers <- ifelse(
    DEGAS.model_type %in% c("ClassClass", "ClassCox"),
    3L,
    0L
  )
  total_layers <- DEGAS.ff_depth + 1L + additional_layers
  thetas <- rlang::list2(
    !!!rlang::set_names(
      lapply(seq_len(total_layers), function(j) py[[glue::glue("Theta{j}")]]),
      glue::glue("theta{j}", j = seq_len(total_layers))
    )
  )

  biases <- rlang::list2(
    !!!rlang::set_names(
      lapply(seq_len(total_layers), function(j) py[[glue::glue("Bias{j}")]]),
      glue::glue("bias{j}", j = seq_len(total_layers))
    )
  )

  methods::new(
    "ccModel",
    Bias = biases,
    Theta = thetas,
    Activation = activation,
    Depth = length(activation),
    Model_type = DEGAS.model_type,
    Architecture = DEGAS.architecture
  )
}
