generateMTLp2 <- function(
  FFdepth = 3L,
  architecture = c("Standard", "DenseNet"),
  ...
) {
  architecture <- match.arg(architecture)
  outlines <- character(0)

  # 单层特殊情况（两种网络类型相同）
  if (FFdepth == 1) {
    return(
      "layerF=add_layer(xs,Fsc,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
    )
  }

  # 标准前馈网络
  if (architecture == "Standard") {
    for (i in seq_len(FFdepth)) {
      if (i == 1) {
        outlines <- c(
          outlines,
          glue::glue(
            "layer1=add_layer(xs,Fsc,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
          )
        )
      } else if (i < FFdepth) {
        outlines <- c(
          outlines,
          glue::glue(
            "layer{i}=add_layer(layer{i-1},hidden_feats,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
          )
        )
      } else {
        # i == FFdepth
        outlines <- c(
          outlines,
          glue::glue(
            "layerF=add_layer(layer{FFdepth-1},hidden_feats,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
          )
        )
      }
    }
    return(paste(outlines, collapse = "\n"))
  }

  #   architecture == "DenseNet"
  for (i in seq_len(FFdepth)) {
    # 1. 构建输入张量连接
    if (i == 1) {
      concat_str <- "xs"
    } else {
      # 动态生成层名称: xs, layer1, layer2, ..., layer{i-1}
      layers_to_concat <- c("xs", stringr::str_c("layer", 1:(i - 1)))
      concat_str <- glue::glue(
        "tf.concat([{stringr::str_c(layers_to_concat, collapse = ', ')}], 1)"
      )
    }

    # 2. 构建输入尺寸字符串
    if (i == 1) {
      inpsz_str <- "Fsc"
    } else {
      inpsz_str <- stringr::str_c(
        "Fsc",
        paste0(rep("+hidden_feats", i - 1), collapse = ""),
        collapse = ""
      )
    }

    # 3. 生成代码行
    if (i < FFdepth) {
      outlines <- c(
        outlines,
        glue::glue(
          "layer{i}=add_layer({concat_str}, {inpsz_str}, hidden_feats, activation_function=tf.sigmoid, dropout_function=True, lambda1=lambda1, keep_prob1=kprob)"
        )
      )
    } else {
      # 最后一层
      outlines <- c(
        outlines,
        glue::glue(
          "layerF=add_layer({concat_str}, {inpsz_str}, hidden_feats, activation_function=tf.sigmoid, dropout_function=True, lambda1=lambda1, keep_prob1=kprob)"
        )
      )
    }
  }

  paste(outlines, collapse = "\n")
}
