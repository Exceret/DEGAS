#' @title ccModel Object
#'
#' @slot Bias A list containing bias terms for each layer of the model.
#' @slot Theta A list of parameter matrices/vectors for the model.
#' @slot Activation A list of activation functions for each layer.
#' @slot Depth Numeric. The number of layers in the model architecture.
#'   Must be a positive integer.
#' @slot Model_type Character. The type of model
#' @slot Architecture Character. Description of the model architecture
#'
#' @examples
#' # Create a simple ccModel object
#' model <- new("ccModel",
#'   Bias = list(c(0.1, 0.2), c(0.3)),
#'   Theta = list(matrix(c(0.5, 0.6), nrow=1), matrix(c(0.7), nrow=1)),
#'   Activation = list("relu", "sigmoid"),
#'   Depth = 2,
#'   Model_type = "neural_network",
#'   Architecture = "feedforward"
#' )
#' @export
#' @family DEGAS_in_SigBridgeR
#' @references Johnson TS, Yu CY, Huang Z, Xu S, Wang T, Dong C, et al. Diagnostic Evidence GAuge of Single cells (DEGAS): a flexible deep transfer learning framework for prioritizing cells in relation to disease. Genome Med. 2022 Feb 1;14(1):11.
#'
setClass(
    "ccModel",
    slots = list(
        Bias = "list",
        Theta = "list",
        Activation = "list",
        Depth = "numeric",
        Model_type = "character",
        Architecture = "character"
    )
)
