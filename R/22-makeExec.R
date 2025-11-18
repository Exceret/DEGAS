# Make python executable for standard (feedforward) implementation

#' @title Generate Python Executable for Neural Network
#' @description
#' Creates a Python executable for feedforward neural network implementation
#' by combining multiple script parts into a complete TensorFlow model.
#'
#' @param tmpDir Temporary directory for script generation
#' @param FFdepth Number of hidden layers in the feedforward network
#' @param model_type Type of neural network model. Must be one of:
#' 'BlankClass', 'ClassBlank', 'BlankCox', 'ClassClass', or 'ClassCox'
#'
#' @return None. Generates Python scripts in the specified temporary directory.
#'
#' @examples
#' \dontrun{
#' # Create a 2-layer classifier model
#' makeExec("/tmp/model/", 2, "ClassClass")
#' }
#'
#' @export
makeExec <- function(tmpDir, FFdepth, model_type) {
    if (
        model_type != 'ClassClass' &&
            model_type != 'ClassCox' &&
            model_type != 'ClassBlank' &&
            model_type != 'BlankClass' &&
            model_type != 'BlankCox'
    ) {
        stop(
            "Please specify either 'BlankClass', 'ClassBlank', 'BlankCox', ClassClass' or 'ClassCox' for the model_type"
        )
    }
    system(paste0('cp ', DEGAS.toolsPath, model_type, 'MTL_p1.py ', tmpDir))
    system(paste0('cp ', DEGAS.toolsPath, model_type, 'MTL_p3.py ', tmpDir))
    outlines = c()
    if (FFdepth == 1) {
        outlines[
            length(outlines) + 1
        ] = "layerF=add_layer(xs,Fsc,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
    } else {
        for (i in 1:FFdepth) {
            if (i == 1) {
                outlines[length(outlines) + 1] = paste0(
                    "layer",
                    as.character(i),
                    "=add_layer(xs,Fsc,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
                )
            } else if (i < FFdepth) {
                outlines[length(outlines) + 1] = paste0(
                    "layer",
                    as.character(i),
                    "=add_layer(layer",
                    as.character(i - 1),
                    ",hidden_feats,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
                )
            } else {
                outlines[length(outlines) + 1] = paste0(
                    "layerF=add_layer(layer",
                    as.character(i - 1),
                    ",hidden_feats,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
                )
            }
        }
    }
    fout = file(paste0(tmpDir, model_type, 'MTL_p2.py'))
    writeLines(outlines, fout)
    close(fout)
    outlines = c()
    outlines[
        length(outlines) + 1
    ] = "#***********************************************************************"
    outlines[length(outlines) + 1] = "# extracting coefficients from TF graph"
    if (model_type == 'ClassClass' || model_type == 'ClassCox') {
        additional_layers = 3
    } else {
        additional_layers = 0
    }
    for (i in 1:(FFdepth + 1 + additional_layers)) {
        if (i == 1) {
            outlines[length(outlines) + 1] = paste0(
                "Theta1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable:0'))[0]"
            )
            outlines[length(outlines) + 1] = paste0(
                "Bias1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_1:0'))[0]"
            )
        } else {
            outlines[length(outlines) + 1] = paste0(
                "Theta",
                as.character(i),
                " = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_",
                as.character((2 * (i - 1))),
                ":0'))[0]"
            )
            outlines[length(outlines) + 1] = paste0(
                "Bias",
                as.character(i),
                " = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_",
                as.character((2 * (i - 1)) + 1),
                ":0'))[0]"
            )
        }
    }
    outlines[
        length(outlines) + 1
    ] = "#***********************************************************************"
    outlines[length(outlines) + 1] = "# Saving model coefficients to files"
    for (i in 1:(FFdepth + 1 + additional_layers)) {
        outlines[length(outlines) + 1] = paste0(
            "np.savetxt(data_folder+'Theta",
            as.character(i),
            ".csv',Theta",
            as.character(i),
            ",delimiter=',')"
        )
        outlines[length(outlines) + 1] = paste0(
            "np.savetxt(data_folder+'Bias",
            as.character(i),
            ".csv',Bias",
            as.character(i),
            ",delimiter=',')"
        )
    }
    outlines[
        length(outlines) + 1
    ] = "with open(data_folder+'Activations.csv','w') as f:"
    for (i in 1:FFdepth) {
        outlines[length(outlines) + 1] = "    f.write('sigmoid\\n')"
    }
    if (model_type == 'BlankCox') {
        outlines[length(outlines) + 1] = "    f.write('sigmoid\\n')"
    } else {
        outlines[length(outlines) + 1] = "    f.write('softmax\\n')"
    }
    if (
        model_type == 'ClassBlank' ||
            model_type == 'BlankClass' ||
            model_type == 'BlankCox'
    ) {
        #outlines[length(outlines)+1] = "    f.write('softmax\\n')"
    } else {
        if (model_type == 'ClassClass') {
            outlines[length(outlines) + 1] = "    f.write('softmax\\n')"
        } else {
            outlines[length(outlines) + 1] = "    f.write('sigmoid\\n')"
        }
        outlines[length(outlines) + 1] = "    f.write('sigmoid\\n')"
        if (model_type == 'ClassClass') {
            outlines[length(outlines) + 1] = "    f.write('softmax\\n')"
        } else {
            outlines[length(outlines) + 1] = "    f.write('sigmoid\\n')"
        }
    }
    fout = file(paste0(tmpDir, model_type, 'MTL_p4.py'))
    writeLines(outlines, fout)
    close(fout)
    system(paste0(
        "cat ",
        tmpDir,
        model_type,
        "MTL_p1.py ",
        tmpDir,
        model_type,
        "MTL_p2.py ",
        tmpDir,
        model_type,
        "MTL_p3.py ",
        tmpDir,
        model_type,
        "MTL_p4.py > ",
        tmpDir,
        model_type,
        "MTL.py"
    ))
}

# Make python executable for densenet implementation

#' @title Generate DenseNet Python Executable
#' @description
#' Creates a Python executable for DenseNet (densely connected neural network)
#' implementation by combining multiple script parts into a complete TensorFlow model.
#'
#' @param tmpDir Temporary directory for script generation
#' @param FFdepth Number of hidden layers in the DenseNet
#' @param model_type Type of neural network model. Must be one of:
#' 'BlankClass', 'ClassBlank', 'BlankCox', 'ClassClass', or 'ClassCox'
#'
#' @return None. Generates Python scripts in the specified temporary directory.
#'
#' @examples
#' \dontrun{
#' # Create a 3-layer DenseNet classifier
#' makeExec2("/tmp/densenet/", 3, "ClassClass")
#' }
#'
#' @export
makeExec2 <- function(tmpDir, FFdepth, model_type) {
    if (
        model_type != 'ClassClass' &&
            model_type != 'ClassCox' &&
            model_type != 'ClassBlank' &&
            model_type != 'BlankClass' &&
            model_type != 'BlankCox'
    ) {
        stop(
            "Please specify either 'BlankClass', 'ClassBlank', 'BlankCox', ClassClass' or 'ClassCox' for the model_type"
        )
    }
    system(paste0('cp ', DEGAS.toolsPath, model_type, 'MTL_p1.py ', tmpDir))
    system(paste0('cp ', DEGAS.toolsPath, model_type, 'MTL_p3.py ', tmpDir))
    outlines = c()
    if (FFdepth == 1) {
        outlines[
            length(outlines) + 1
        ] = "layerF=add_layer(xs,Fsc,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
    } else {
        inpsz_str = "Fsc"
        for (i in 1:FFdepth) {
            if (i == 1) {
                outlines[length(outlines) + 1] = paste0(
                    "layer",
                    as.character(i),
                    "=add_layer(xs,",
                    inpsz_str,
                    ",hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
                )
            } else if (i < FFdepth) {
                inpsz_str = paste0(inpsz_str, "+hidden_feats")
                layer_str = "tf.concat([xs"
                for (j in 1:i) {
                    if (j == i) {
                        layer_str = paste0(layer_str, '],1)')
                    } else {
                        layer_str = paste0(layer_str, ',layer', as.character(j))
                    }
                }
                outlines[length(outlines) + 1] = paste0(
                    "layer",
                    as.character(i),
                    "=add_layer(",
                    layer_str,
                    ",",
                    inpsz_str,
                    ",hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
                )
            } else {
                inpsz_str = paste0(inpsz_str, "+hidden_feats")
                layer_str = "tf.concat([xs"
                for (j in 1:FFdepth) {
                    if (j == i) {
                        layer_str = paste0(layer_str, '],1)')
                    } else {
                        layer_str = paste0(layer_str, ',layer', as.character(j))
                    }
                }
                outlines[length(outlines) + 1] = paste0(
                    "layerF=add_layer(",
                    layer_str,
                    ",",
                    inpsz_str,
                    ",hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
                )
            }
        }
    }
    fout = file(paste0(tmpDir, model_type, 'MTL_p2.py'))
    writeLines(outlines, fout)
    close(fout)
    outlines = c()
    outlines[
        length(outlines) + 1
    ] = "#***********************************************************************"
    outlines[length(outlines) + 1] = "# extracting coefficients from TF graph"
    if (model_type == 'ClassClass' || model_type == 'ClassCox') {
        additional_layers = 3
    } else {
        additional_layers = 0
    }
    for (i in 1:(FFdepth + 1 + additional_layers)) {
        if (i == 1) {
            outlines[length(outlines) + 1] = paste0(
                "Theta1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable:0'))[0]"
            )
            outlines[length(outlines) + 1] = paste0(
                "Bias1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_1:0'))[0]"
            )
        } else {
            outlines[length(outlines) + 1] = paste0(
                "Theta",
                as.character(i),
                " = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_",
                as.character((2 * (i - 1))),
                ":0'))[0]"
            )
            outlines[length(outlines) + 1] = paste0(
                "Bias",
                as.character(i),
                " = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_",
                as.character((2 * (i - 1)) + 1),
                ":0'))[0]"
            )
        }
    }
    outlines[
        length(outlines) + 1
    ] = "#***********************************************************************"
    outlines[length(outlines) + 1] = "# Saving model coefficients to files"
    for (i in 1:(FFdepth + 1 + additional_layers)) {
        outlines[length(outlines) + 1] = paste0(
            "np.savetxt(data_folder+'Theta",
            as.character(i),
            ".csv',Theta",
            as.character(i),
            ",delimiter=',')"
        )
        outlines[length(outlines) + 1] = paste0(
            "np.savetxt(data_folder+'Bias",
            as.character(i),
            ".csv',Bias",
            as.character(i),
            ",delimiter=',')"
        )
    }
    outlines[
        length(outlines) + 1
    ] = "with open(data_folder+'Activations.csv','w') as f:"
    for (i in 1:FFdepth) {
        outlines[length(outlines) + 1] = "    f.write('sigmoid\\n')"
    }
    if (model_type == 'BlankCox') {
        outlines[length(outlines) + 1] = "    f.write('sigmoid\\n')"
    } else {
        outlines[length(outlines) + 1] = "    f.write('softmax\\n')"
    }
    if (
        model_type == 'ClassBlank' ||
            model_type == 'BlankClass' ||
            model_type == 'BlankCox'
    ) {
        #outlines[length(outlines)+1] = "    f.write('softmax\\n')"
    } else {
        if (model_type == 'ClassClass') {
            outlines[length(outlines) + 1] = "    f.write('softmax\\n')"
        } else {
            outlines[length(outlines) + 1] = "    f.write('sigmoid\\n')"
        }
        outlines[length(outlines) + 1] = "    f.write('sigmoid\\n')"
        if (model_type == 'ClassClass') {
            outlines[length(outlines) + 1] = "    f.write('softmax\\n')"
        } else {
            outlines[length(outlines) + 1] = "    f.write('sigmoid\\n')"
        }
    }
    fout = file(paste0(tmpDir, model_type, 'MTL_p4.py'))
    writeLines(outlines, fout)
    close(fout)
    system(paste0(
        "cat ",
        tmpDir,
        model_type,
        "MTL_p1.py ",
        tmpDir,
        model_type,
        "MTL_p2.py ",
        tmpDir,
        model_type,
        "MTL_p3.py ",
        tmpDir,
        model_type,
        "MTL_p4.py > ",
        tmpDir,
        model_type,
        "MTL.py"
    ))
}
