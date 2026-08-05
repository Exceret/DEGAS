"""
Modern DEGAS implementation for reticulate.

Main API:
    run_degas(...)
    train_degas(...)
    predict_degas(...)
    DEGASTensorFlow

Compatibility-like names:
    preprocessCounts
    toCorrCoeff
    knnSmooth
    runCCMTLBag
    predClassBag
"""

from .api import (
    DEGASTensorFlow,
    run_degas,
    train_degas,
    predict_degas,
    runCCMTLBag,
    predClassBag,
)

from .preprocessing import (
    normFunc,
    scaleFunc,
    preprocessCounts,
    align_genes,
    scale_expression,
)

from .postprocessing import (
    centerFunc,
    toCorrCoeff,
    knnSmooth,
)

__all__ = [
    "DEGASTensorFlow",
    "run_degas",
    "train_degas",
    "predict_degas",
    "runCCMTLBag",
    "predClassBag",
    "normFunc",
    "scaleFunc",
    "preprocessCounts",
    "align_genes",
    "scale_expression",
    "centerFunc",
    "toCorrCoeff",
    "knnSmooth",
]
