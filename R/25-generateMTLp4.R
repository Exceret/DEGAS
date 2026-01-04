generateMTLp4 <- function(
  FFdepth = 3L,
  model_type = c(
    'ClassClass',
    'ClassCox',
    'ClassBlank',
    'BlankClass',
    'BlankCox'
  ),
  ...
) {
  model_type <- match.arg(model_type)
  outlines <- character(0)

  separator <- "# ***********************************************************************"

  # 确定额外层数
  additional_layers <- ifelse(model_type %in% c('ClassClass', 'ClassCox'), 3, 0)
  total_layers <- FFdepth + 1 + additional_layers

  # 系数提取部分
  outlines <- c(outlines, separator, "# extracting coefficients from TF graph")

  for (i in seq_len(total_layers)) {
    if (i == 1) {
      outlines <- c(
        outlines,
        "Theta1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, 'Variable:0'))[0]",
        "Bias1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, 'Variable_1:0'))[0]"
      )
    } else {
      var_idx_w <- 2 * (i - 1)
      var_idx_b <- var_idx_w + 1
      outlines <- c(
        outlines,
        glue::glue(
          "Theta{i} = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, 'Variable_{var_idx_w}:0'))[0]"
        ),
        glue::glue(
          "Bias{i} = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, 'Variable_{var_idx_b}:0'))[0]"
        )
      )
    }
  }

  #   # 保存系数部分
  #   outlines <- c(outlines, separator, "# Saving model coefficients to files")

  #   for (i in seq_len(total_layers)) {
  #     outlines <- c(
  #       outlines,
  #       glue::glue(
  #         "np.savetxt(data_folder+ 'Theta{i}.csv', Theta{i}, delimiter=',')"
  #       ),
  #       glue::glue(
  #         "np.savetxt(data_folder+ 'Bias{i}.csv', Bias{i}, delimiter=',')"
  #       )
  #     )
  #   }

  # 激活函数记录部分
  activation_str <- "\nactivation: list = ["

  # 隐藏层激活函数
  for (i in seq_len(FFdepth)) {
    if (i == 1) {
      activation_str <- paste(activation_str, "'sigmoid'", sep = "")
    } else {
      activation_str <- paste(activation_str, "'sigmoid'", sep = ", ")
    }
  }

  # 输出层激活函数
  if (model_type == 'BlankCox') {
    activation_str <- paste(activation_str, "'sigmoid'", sep = ", ")
  } else {
    activation_str <- paste(activation_str, "'softmax'", sep = ", ")
  }

  # 额外输出层 (根据模型类型)
  if (!(model_type %in% c('ClassBlank', 'BlankClass', 'BlankCox'))) {
    if (model_type == 'ClassClass') {
      activation_str <- paste(activation_str, "'softmax'", sep = ", ")
    } else {
      activation_str <- paste(activation_str, "'sigmoid'", sep = ", ")
    }
    activation_str <- paste(activation_str, "'sigmoid'", sep = ", ")

    if (model_type == 'ClassClass') {
      activation_str <- paste(activation_str, "'softmax'", sep = ", ")
    } else {
      activation_str <- paste(activation_str, "'sigmoid'", sep = ", ")
    }
  }
  activation_str <- paste0(activation_str, "]")
  outlines <- c(outlines, activation_str)

  paste(outlines, collapse = "\n")
}
