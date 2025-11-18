#***************************************************************
#***************************************************************
# Single Cell Transfer Learning Toolbox
# Travis S Johnson and Zhi Huang
# Compatible with OSX and Linux

#***************************************************************
#***************************************************************
# Initializing the package

#' @title Initialize DEGAS
#'
#' @description Initialize the DEGAS package
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' initDEGAS()
#' }
#'
#' @export
initDEGAS <- function() {
    DEGAS.pyloc <<- "python3"
    DEGAS.toolsPath <<- paste0(.libPaths()[1], "/DEGAS/tools/")
    DEGAS.train_steps <<- 2000
    DEGAS.scbatch_sz <<- 200
    DEGAS.patbatch_sz <<- 50
    DEGAS.hidden_feats <<- 50
    DEGAS.do_prc <<- 0.5
    DEGAS.lambda1 <<- 3.0
    DEGAS.lambda2 <<- 3.0
    DEGAS.lambda3 <<- 3.0
    DEGAS.seed <<- "NULL"
}

#***************************************************************
#***************************************************************
# General utility functions

# Return filename extension
#' @keywords internal
getExtension <- function(file) {
    ex <- strsplit(basename(file), split = "\\.")[[1]]
    return(ex[length(ex)])
}

# Manually reset the python path
#' @keywords internal
setPython <- function(path2python) {
    DEGAS.pyloc <<- path2python
    #Sys.setenv(PATH = paste(c(path2python,Sys.getenv("PATH")),collapse = .Platform$path.sep))
}

# Manually reset the number of training steps
#' @keywords internal
set_training_steps <- function(inp) {
    DEGAS.train_steps <<- inp
}

# Manually reset the single cell batch size
#' @keywords internal
set_single_cell_batch_size <- function(inp) {
    DEGAS.scbatch_sz <<- inp
}

# Manually reset the single patient batch size
#' @keywords internal
set_patient_batch_size <- function(inp) {
    DEGAS.patbatch_sz <<- inp
}

# Manually reset the number of hidden features
#' @keywords internal
set_hidden_feature_number <- function(inp) {
    DEGAS.hidden_feats <<- inp
}

# Manually reset the dropout keep percentage (the percentage of nodes to keep)
#' @keywords internal
set_dropout_keep_fraction <- function(inp) {
    DEGAS.do_prc <<- inp
}

# Manually reset the L2 regularization term (lambda 3)
#' @keywords internal
set_l2_regularization_term <- function(inp) {
    DEGAS.lambda1 <<- inp
}

# Manually reset the patient loss term (lambda 1)
#' @keywords internal
set_patient_loss_term <- function(inp) {
    DEGAS.lambda2 <<- inp
}

# Manually reset the MMD loss term (lambda 2)
#' @keywords internal
set_MMD_loss_term <- function(inp) {
    DEGAS.lambda3 <<- inp
}

# Manually reset the seed
#' @title Set seed
#' @param inp seed
#' @export
set_seed_term <- function(inp) {
    if (is.null(inp)) {
        DEGAS.seed <<- "NULL"
    } else if (is.numeric(inp) & length(inp) == 1) {
        DEGAS.seed <<- floor(inp)
    } else {
        message("ERROR: Input proper seed (NULL or integer)")
        message("Non-integer values will equal floor(value)")
    }
}

#***************************************************************
#***************************************************************a

# zscore normalization
#' @keywords internal
normFunc <- function(x) {
    return((x - Matrix::mean(x, na.rm = T)) / (stats::sd(x, na.rm = T) + 1e-3))
}

# scaling from 0-1
#' @keywords internal
scaleFunc <- function(x) {
    return(
        (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T) + 1e-3)
    )
}

# Preprocess count data
#' @keywords internal
normalizeScale <- function(X) {
    return(Matrix::t(apply(
        Matrix::t(apply(Matrix::as.matrix(Matrix::t(X)), 1, normFunc)),
        1,
        scaleFunc
    )))
}

#' @keywords internal
preprocessCounts <- function(X) {
    return(normalizeScale(1.5^log2(X + 1)))
}

#' @keywords internal
preprocessAllCounts <- function(sc.dat, pt.dat) {
    ks.d = c(0, 0, 0, 0)
    names(ks.d) = c("nonenone", "log2none", "nonelog2", "log2log2")
    tmp = suppressWarnings(stats::ks.test(
        preprocessCounts(sc.dat),
        preprocessCounts(pt.dat)
    ))
    ks.d["nonenone"] = tmp$statistic
    tmp = suppressWarnings(stats::ks.test(
        preprocessCounts(log2(sc.dat + 1)),
        preprocessCounts(pt.dat)
    ))
    ks.d["log2none"] = tmp$statistic
    tmp = suppressWarnings(stats::ks.test(
        preprocessCounts(sc.dat),
        preprocessCounts(log2(pt.dat + 1))
    ))
    ks.d["nonelog2"] = tmp$statistic
    tmp = suppressWarnings(stats::ks.test(
        preprocessCounts(log2(sc.dat + 1)),
        preprocessCounts(log2(pt.dat + 1))
    ))
    ks.d["log2log2"] = tmp$statistic
    message(names(ks.d)[ks.d == min(ks.d)])
    if (names(ks.d)[ks.d == min(ks.d)] == "nonenone") {
        return(list(
            scDat = preprocessCounts(sc.dat),
            patDat = preprocessCounts(pt.dat)
        ))
    } else if (names(ks.d)[ks.d == min(ks.d)] == "log2none") {
        return(list(
            scDat = preprocessCounts(log2(sc.dat + 1)),
            patDat = preprocessCounts(pt.dat)
        ))
    } else if (names(ks.d)[ks.d == min(ks.d)] == "nonelog2") {
        return(list(
            scDat = preprocessCounts(sc.dat),
            patDat = preprocessCounts(log2(pt.dat + 1))
        ))
    } else {
        return(list(
            scDat = preprocessCounts(log2(sc.dat + 1)),
            patDat = preprocessCounts(log2(pt.dat + 1))
        ))
    }
}

# center to 0
#' @keywords internal
centerFunc <- function(x) {
    return(x - Matrix::mean(x, na.rm = T))
}

# Activation functions and utilities

# Sigmoid activation function
#' @keywords internal
sigmoid <- function(x) {
    1 / (1 + exp(-x))
}

# Log sum exp transformation (for softmax)
#' @keywords internal
logsumexp <- function(x) {
    y = max(x)
    y + log(sum(exp(x - y)))
}

# Softmax activation function
#' @keywords internal
softmax <- function(X) {
    return(Matrix::t(apply(X, 1, function(x) exp(x - logsumexp(x)))))
}


# List of label names to a onehot matrix with labels as column names
#' @keywords internal
toOneHot <- function(labels) {
    labs = unique(labels)
    out = matrix(0, length(labels), length(labs))
    colnames(out) = labs
    row.names(out) = row.names(labels)
    for (i in 1:length(labels)) {
        out[i, labels[i]] = 1
    }
    return(out)
}

# Convert matrix of output weights to max value for each row (row max = 1 and not row max = 0)
#' @keywords internal
probtoOneHot <- function(probMat) {
    idx = apply(probMat, 1, function(x) Matrix::which(x == max(x)))
    probMat = probMat * 0
    for (i in 1:length(idx)) {
        probMat[i, idx[i]] = 1
    }
    return(probMat)
}

# Predict patient class from proportions of single cell classes
#' @keywords internal
predPatClassFromSCClass <- function(ccModel1, Exp) {
    Z1 = sigmoid(sweep((Exp %*% ccModel4@Theta4), 2, ccModel1@Bias4, '+'))
    return(softmax(sweep((Z1 %*% ccModel1@Theta5), 2, ccModel1@Bias5, '+')))
}

#***************************************************************
# Other functions

# Generates sets for k fold cross validation
#' @keywords internal
splitKfoldCV <- function(N, k) {
    if (k < 3) {
        stop("Please use 3 or more folds")
    }
    Idx = as.numeric(sample(1:N, N, replace = FALSE))
    sz = rep(floor(N / k), k)
    rem = N - sum(sz)
    if (rem > 0) {
        cntr = 0
        for (i in 1:rem) {
            if (cntr == k) {
                cntr = 1
            } else {
                cntr = cntr + 1
            }
            sz[cntr] = sz[cntr] + 1
        }
    }
    cntr = 0
    grpIdx = list()
    for (i in 1:k) {
        grpIdx[[i]] = Idx[(cntr + 1):(cntr + sz[i])]
        cntr = cntr + sz[i]
    }
    return(grpIdx)
}

# Get a feature vector from a dataframe
#' @keywords internal
getFeat = function(vec, df, colm, colo) {
    tmp = vec
    for (i in 1:length(vec)) {
        tmp[i] = df[Matrix::which(df[, colm] == vec[i])[1], colo]
    }
    return(tmp)
}

# returns duplicate row names to remove
#' @keywords internal
remDupIdx <- function(X, dup_rnames, rnames) {
    rem = c()
    for (dup_rname in dup_rnames) {
        tmp = Matrix::which(rnames == dup_rname)
        Xtmp = X[tmp, ]
        Xmean = Matrix::rowMeans(Matrix::as.matrix(Xtmp[, 2:dim(Xtmp)[2]]))
        Xmean[is.na(Xmean)] = 0
        Xmean = abs(Xmean)
        rem = c(rem, tmp[Matrix::which(Xmean != max(Xmean, na.rm = TRUE))])
    }
    return(rem)
}

#***************************************************************
# Post-processing functions

# Quantile normalization
#' @keywords internal
quantNorm <- function(
    df,
    center = 'median',
    rescale = TRUE,
    rescale_mult = 1e4
) {
    df_rank <- apply(df, 2, rank, ties.method = "min")
    df_sorted <- data.frame(apply(df, 2, sort))
    df_mean <- apply(df_sorted, 1, mean)
    index_to_mean <- function(my_index, my_mean) {
        return(my_mean[my_index])
    }
    df_final <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)
    rownames(df_final) <- rownames(df)
    meds = eval(parse(text = paste0("apply(df_final, 2, ", center, ")")))
    df_final = sweep(df_final, 2, meds, '-')
    if (rescale) {
        df_final[df_final > 0] = eval(parse(
            text = paste0("log2(", rescale_mult, "*df_final[df_final>0])")
        ))
        df_final[df_final < 0] = eval(parse(
            text = paste0("-log2(-", rescale_mult, "*df_final[df_final<0])")
        ))
        meds = eval(parse(text = paste0("apply(df_final, 2, ", center, ")")))
        df_final = sweep(df_final, 2, meds, '-')
    }
    return(df_final)
}

# Return euclidean distance between two points
#' @keywords internal
euclDist <- function(loc1, loc2) {
    return(sqrt(sum((loc1 - loc2)^2)))
}

# Return a matrix of all pairwise distances
#' @keywords internal
pairDist <- function(locs) {
    N = dim(locs)[1]
    out = matrix(NA, N, N)
    for (i in 1:N) {
        for (j in 1:i) {
            out[i, j] = out[j, i] = euclDist(locs[i, ], locs[j, ])
        }
    }
    return(out)
}

# Return k-nearest-neighbor smoothed probabilites
#' @keywords internal
knnSmooth <- function(probs, locs, k = 5) {
    out = probs
    dists = pairDist(locs)
    if (class(probs)[1] == "numeric") {
        N = length(probs)
        for (i in 1:N) {
            idx = order(dists[i, ])
            out[i] = Matrix::mean(probs[idx[1:k]], na.rm = TRUE)
        }
    } else {
        N = dim(locs)[1]
        for (i in 1:N) {
            idx = order(dists[i, ])
            out[i, ] = Matrix::colMeans(probs[idx[1:k], ], na.rm = TRUE)
        }
    }
    return(out)
}

# return random sample (of s) which is evenly distributed across sample groups (g)
# where each group has n samples.
# Note: If a group has less than n samples, then all samples in that group are used.
#' @keywords internal
evenSamp <- function(s, g, n) {
    groups = unique(g)
    out = list()
    for (group in groups) {
        if (sum(g == group) >= n) {
            out[[group]] = sample(s[g == group], n, replace = FALSE)
        } else if (sum(g == group) > 0) {
            out[[group]] = sample(
                s[g == group],
                sum(g == group),
                replace = FALSE
            )
        } else {
            #Adding nothing
        }
    }
    out = unlist(out)
    return(out)
}

# Convert DEGAS output [0,1] to an association [-1,1]
#' @keywords internal
toCorrCoeff <- function(probs) {
    k = dim(probs)[2]
    if (k < 2 || is.null(k)) {
        k = 2
    }
    l = 2
    return(2 * ((probs - 1 / k) / (l - l / k) + 1 / l) - 1)
}

#***************************************************************
# Functions for atlas level datasets and additional bootstrapping statistics

# Running DEGAS for large, atlas level, datasets
# ! Future update
#' @keywords internal
runDEGASatlas <- function(
    scDat,
    scLab,
    patDat,
    patLab,
    tmpDir,
    model_type,
    architecture,
    FFdepth,
    Bagdepth,
    Nsubsample,
    seed
) {
    folds = splitKfoldCV(dim(scDat)[1], floor(dim(scDat)[1] / Nsubsample))
    ccModel_out = list()
    for (f in folds) {
        initDEGAS()
        set_seed_term(seed)
        ccModel_tmp = runCCMTLBag(
            scDat[f, ],
            scLab[f, ],
            patDat,
            patLab,
            tmpDir,
            model_type,
            architecture,
            FFdepth,
            Bagdepth
        )
        ccModel_out = c(ccModel_out, ccModel_tmp)
    }
    return(ccModel_out)
}

# Predict quantiles and statistics for Impressions
#' @keywords internal
getCI <- function(x) {
    tmp = c(Matrix::mean(x), stats::quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)))
    names(tmp) = c("Mean", "Prc5", "Prc25", "Prc50", "Prc75", "Prc95")
    tmp["IQR"] = tmp["Prc75"] - tmp["Prc25"]
    tmp["CI95diff"] = tmp["Prc95"] - tmp["Prc5"]
    return(tmp)
}

#' @keywords internal
predClassBagCI <- function(ccModel, Exp, scORpat) {
    out = list()
    for (i in 1:length(ccModel)) {
        out[[i]] <- predClass(ccModel[[i]], Exp, scORpat)
    }
    n_reps = length(out)
    n_labs = dim(out[[1]])[2]
    n_cells = dim(out[[1]])[1]
    out = do.call(cbind, out)
    tmp = list()
    for (i in 1:n_labs) {
        message(i)
        tmp[[i]] = Matrix::t(apply(
            out[, seq(i, dim(out)[2], n_labs)],
            1,
            function(x) {
                getCI(x)
            }
        ))
        colnames(tmp[[i]]) = paste0(colnames(tmp[[1]]), "_", as.character(i))
    }
    out = do.call(cbind, tmp)
    return(out)
}

#  knnSmoothing for large, atlas level, datasets
#' @keywords internal
knnSmoothAtlas <- function(sc_seurat, preds, k, n_split) {
    folds = splitKfoldCV(dim(sc_seurat)[2], floor(dim(sc_seurat)[2] / n_split))
    cor_list = list()
    i = 0
    for (f in folds) {
        i = i + 1
        samp = f
        umap_coords = sc_seurat@reductions$umap@cell.embeddings[samp, ]
        predsSmoothed = knnSmooth(
            Matrix::as.matrix(preds[samp, ]),
            Matrix::as.matrix(umap_coords),
            k
        )
        RespCor = toCorrCoeff(scaleFunc(predsSmoothed))
        df = data.frame(
            Idx = f,
            UMAP_1 = umap_coords[, 1],
            UMAP_2 = umap_coords[, 2],
            RespCor = RespCor,
            Cluster = sc_seurat$seurat_clusters[samp]
        )
        cor_list[[i]] = df
    }
    cor_list_df = do.call(rbind, cor_list)
    return(cor_list_df[rownames(preds), ])
}
