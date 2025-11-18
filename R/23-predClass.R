# Make predictions based on a trained DEGAS model

#' @title Make Predictions Using Trained DEGAS Model
#' @description
#' Makes class predictions based on a trained DEGAS model. Automatically
#' routes to the appropriate prediction function based on the model architecture.
#'
#' @param ccModel1 A trained DEGAS model object
#' @param Exp Expression matrix for prediction
#' @param scORpat Indicator for single-cell or patient data
#'
#' @return Model predictions
#'
#' @examples
#' \dontrun{
#' # Make predictions using a trained model
#' predictions <- predClass(trained_model, expression_matrix, "sc")
#' }
#'
#' @export
predClass <- function(ccModel1, Exp, scORpat) {
    if (ccModel1@Architecture == "DenseNet") {
        return(predClass2(ccModel1, Exp, scORpat))
    } else if (ccModel1@Architecture == "Standard") {
        return(predClass1(ccModel1, Exp, scORpat))
    } else {
        stop("Incorrect architecture argument")
    }
}

# Prediction from trained standard architecture model

#' @title Prediction from Standard Architecture Model
#' @description
#' Makes predictions using a trained DEGAS model with standard feedforward architecture.
#'
#' @param ccModel1 A trained DEGAS model object
#' @param Exp Expression matrix for prediction
#' @param scORpat Indicator for single-cell ("sc") or patient ("pat") data
#'
#' @return Model predictions
#'
#' @export
predClass1 <- function(ccModel1, Exp, scORpat) {
    Z = Exp
    rm(Exp)
    if (
        ccModel1@Model_type == 'BlankClass' ||
            ccModel1@Model_type == 'ClassBlank' ||
            ccModel1@Model_type == 'ClassBlank'
    ) {
        for (i in 1:(ccModel1@Depth)) {
            calcZ = paste0(
                ccModel1@Activation[[i]],
                "(sweep((as.matrix(Z) %*% ccModel1@Theta[[",
                as.character(i),
                "]]),2,ccModel1@Bias[[",
                as.character(i),
                "]],'+'))"
            )
            Z = eval(parse(text = calcZ))
        }
        return(Z)
    } else {
        for (i in 1:(ccModel1@Depth - 4)) {
            calcZ = paste0(
                ccModel1@Activation[[i]],
                "(sweep((as.matrix(Z) %*% ccModel1@Theta[[",
                as.character(i),
                "]]),2,ccModel1@Bias[[",
                as.character(i),
                "]],'+'))"
            )
            Z = eval(parse(text = calcZ))
        }
    }
    if (toupper(scORpat) == 'SC') {
        calcPred = paste0(
            ccModel1@Activation[[ccModel1@Depth - 3]],
            "(sweep((Z %*% ccModel1@Theta[[",
            as.character(ccModel1@Depth - 3),
            "]]),2,ccModel1@Bias[[",
            as.character(ccModel1@Depth - 3),
            "]],'+'))"
        )
    } else if (toupper(scORpat) == 'PAT') {
        calcPred = paste0(
            ccModel1@Activation[[ccModel1@Depth - 2]],
            "(sweep((Z %*% ccModel1@Theta[[",
            as.character(ccModel1@Depth - 2),
            "]]),2,ccModel1@Bias[[",
            as.character(ccModel1@Depth - 2),
            "]],'+'))"
        )
    } else {
        stop("Incorrect prediction argument. Please use 'sc' or 'pat'")
    }
    return(eval(parse(text = calcPred)))
}

#' @title Prediction from DenseNet Architecture Model
#' @description
#' Makes predictions using a trained DEGAS model with DenseNet architecture.
#'
#' @param ccModel1 A trained DEGAS model object
#' @param Exp Expression matrix for prediction
#' @param scORpat Indicator for single-cell ("sc") or patient ("pat") data
#'
#' @return Model predictions
#'
#' @export
predClass2 <- function(ccModel1, Exp, scORpat) {
    Z = Exp
    rm(Exp)
    if (
        ccModel1@Model_type == 'BlankClass' ||
            ccModel1@Model_type == 'ClassBlank' ||
            ccModel1@Model_type == 'BlankCox'
    ) {
        for (i in 1:(ccModel1@Depth)) {
            calcZ = paste0(
                ccModel1@Activation[[i]],
                "(sweep((as.matrix(Z) %*% ccModel1@Theta[[",
                as.character(i),
                "]]),2,ccModel1@Bias[[",
                as.character(i),
                "]],'+'))"
            )
            if (i < ccModel1@Depth - 1) {
                Z = cbind(Z, eval(parse(text = calcZ)))
            } else {
                Z = eval(parse(text = calcZ))
            }
        }
        return(Z)
    } else {
        for (i in 1:(ccModel1@Depth - 4)) {
            calcZ = paste0(
                ccModel1@Activation[[i]],
                "(sweep((as.matrix(Z) %*% ccModel1@Theta[[",
                as.character(i),
                "]]),2,ccModel1@Bias[[",
                as.character(i),
                "]],'+'))"
            )
            if (i < ccModel1@Depth - 4) {
                Z = cbind(Z, eval(parse(text = calcZ)))
            } else {
                Z = eval(parse(text = calcZ))
            }
        }
    }
    if (toupper(scORpat) == 'SC') {
        calcPred = paste0(
            ccModel1@Activation[[ccModel1@Depth - 3]],
            "(sweep((Z %*% ccModel1@Theta[[",
            as.character(ccModel1@Depth - 3),
            "]]),2,ccModel1@Bias[[",
            as.character(ccModel1@Depth - 3),
            "]],'+'))"
        )
    } else if (toupper(scORpat) == 'PAT') {
        calcPred = paste0(
            ccModel1@Activation[[ccModel1@Depth - 2]],
            "(sweep((Z %*% ccModel1@Theta[[",
            as.character(ccModel1@Depth - 2),
            "]]),2,ccModel1@Bias[[",
            as.character(ccModel1@Depth - 2),
            "]],'+'))"
        )
    } else {
        stop("Incorrect prediction argument. Please use 'sc' or 'pat'")
    }
    return(eval(parse(text = calcPred)))
}
