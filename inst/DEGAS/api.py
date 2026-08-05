import os

os.environ.setdefault("KERAS_BACKEND", "tensorflow")

from collections.abc import Mapping

import numpy as np
import pandas as pd
import tensorflow as tf
import keras

from .model import DenseNetDEGAS
from .losses import cox_ph_loss, mmd_rbf_loss
from .preprocessing import align_genes, scale_expression
from .feature_selection import DEGASFeatureSelector
from .ts_print import ts_print


class DEGASTensorFlow:
    """
    Modern TensorFlow/Keras DEGAS implementation.

    Input convention
    ----------------
    sc_matrix:
        cells × genes

    bulk_matrix:
        bulk samples × genes

    Cox bulk_labels:
        pandas.DataFrame with columns:
            time
            status

    Classification bulk_labels:
        dict:
            {bulk_sample_name: 0/1}
        or array-like labels.
    """

    def __init__(
        self,
        patient_task="auto",
        n_bootstrap=5,
        training_steps=2000,
        cell_batch_size=200,
        bulk_batch_size=50,
        hidden_units=50,
        n_layers=3,
        dropout_rate=0.5,
        lambda_cell=2.0,
        lambda_patient=3.0,
        lambda_mmd=3.0,
        lambda_l2=3.0,
        l2_base=1e-4,
        learning_rate=1e-3,
        top_var_frac=0.2,
        max_features=1000,
        min_features=50,
        random_state=123,
        verbose=True,
        log_transform=False,
        cox_l1_ratio=0.5,
        cox_n_alphas=100,
        cox_alpha_min_ratio=0.01,
        cox_max_iter=100000,
        cox_tol=1e-7,
        cox_top_n=None,
    ):
        self.patient_task = patient_task

        self.n_bootstrap = int(n_bootstrap)
        self.training_steps = int(training_steps)
        self.cell_batch_size = int(cell_batch_size)
        self.bulk_batch_size = int(bulk_batch_size)

        self.hidden_units = int(hidden_units)
        self.n_layers = int(n_layers)
        self.dropout_rate = float(dropout_rate)

        self.lambda_cell = float(lambda_cell)
        self.lambda_patient = float(lambda_patient)
        self.lambda_mmd = float(lambda_mmd)
        self.lambda_l2 = float(lambda_l2)
        self.l2_base = float(l2_base)
        self.learning_rate = float(learning_rate)

        self.top_var_frac = float(top_var_frac)
        self.max_features = None if max_features is None else int(max_features)
        self.min_features = int(min_features)

        self.random_state = int(random_state)
        self.verbose = bool(verbose)
        self.log_transform = bool(log_transform)

        self.cox_l1_ratio = float(cox_l1_ratio)
        self.cox_n_alphas = int(cox_n_alphas)
        self.cox_alpha_min_ratio = float(cox_alpha_min_ratio)
        self.cox_max_iter = int(cox_max_iter)
        self.cox_tol = float(cox_tol)
        self.cox_top_n = None if cox_top_n is None else int(cox_top_n)

        self.models_ = []
        self.selector_ = None

        self.patient_task_ = None
        self.common_genes_ = None
        self.selected_genes_ = None
        self.selected_gene_indices_ = None

        self.cell_names_ = None
        self.bulk_sample_names_ = None

        self.patient_classes_ = None
        self.patient_class_to_index_ = None
        self.disease_class_index_ = None

        self.cell_classes_ = None
        self.cell_class_to_index_ = None

    # ------------------------------------------------------------------
    # Label preparation
    # ------------------------------------------------------------------

    def _infer_patient_task(self, bulk_labels):
        if self.patient_task != "auto":
            return self.patient_task

        if isinstance(bulk_labels, pd.DataFrame):
            cols = [str(c).lower() for c in bulk_labels.columns]
            if "time" in cols and "status" in cols:
                return "cox"

        if isinstance(bulk_labels, Mapping):
            keys = set([str(k).lower() for k in bulk_labels.keys()])
            if "time" in keys and "status" in keys:
                return "cox"
            return "classification"

        return "classification"

    def _prepare_bulk_labels(self, bulk_labels, bulk_sample_names, n_bulk):
        task = self._infer_patient_task(bulk_labels)

        if task == "cox":
            if isinstance(bulk_labels, pd.DataFrame):
                df = bulk_labels.copy()
                colmap = {str(c).lower(): c for c in df.columns}

                if "time" not in colmap or "status" not in colmap:
                    raise ValueError(
                        "Cox label DataFrame must contain time and status."
                    )

                if bulk_sample_names is not None:
                    bulk_sample_names = list(map(str, bulk_sample_names))
                    if set(bulk_sample_names).issubset(set(map(str, df.index))):
                        df.index = list(map(str, df.index))
                        df = df.loc[bulk_sample_names]
                    elif len(df) != n_bulk:
                        raise ValueError(
                            "bulk_sample_names not found in Cox label index, "
                            "and label length != number of bulk samples."
                        )

                if len(df) != n_bulk:
                    raise ValueError(
                        "Cox label rows must equal number of bulk samples."
                    )

                time = df[colmap["time"]].values.astype(np.float32)
                status = df[colmap["status"]].values.astype(np.float32)

            elif (
                isinstance(bulk_labels, Mapping)
                and "time" in bulk_labels
                and "status" in bulk_labels
            ):
                time = np.asarray(bulk_labels["time"], dtype=np.float32)
                status = np.asarray(bulk_labels["status"], dtype=np.float32)

                if len(time) != n_bulk or len(status) != n_bulk:
                    raise ValueError("Cox time/status length must equal bulk samples.")

            else:
                raise ValueError(
                    "For Cox task, bulk_labels must be pandas.DataFrame "
                    "with time/status or dict with time/status."
                )

            if np.any(time <= 0):
                raise ValueError("Cox survival time must be positive.")

            if not np.all(np.isin(status, [0, 1])):
                raise ValueError("Cox status must contain 0/1.")

            return "cox", {"time": time, "status": status}, 1

        elif task == "classification":
            if isinstance(bulk_labels, Mapping):
                if bulk_sample_names is None:
                    raise ValueError(
                        "For dict classification labels, bulk_sample_names is required."
                    )

                bulk_sample_names = list(map(str, bulk_sample_names))

                y_raw = []
                for s in bulk_sample_names:
                    if s in bulk_labels:
                        y_raw.append(bulk_labels[s])
                    elif str(s) in bulk_labels:
                        y_raw.append(bulk_labels[str(s)])
                    else:
                        raise KeyError(f"Sample {s} not found in bulk_labels.")
                y_raw = np.asarray(y_raw)

            else:
                y_raw = np.asarray(bulk_labels)

                if y_raw.ndim == 2:
                    y_raw = np.argmax(y_raw, axis=1)

                if len(y_raw) != n_bulk:
                    raise ValueError("Classification label length != bulk samples.")

            classes = sorted(pd.unique(y_raw))
            class_to_index = {c: i for i, c in enumerate(classes)}
            y = np.asarray([class_to_index[v] for v in y_raw], dtype=np.int32)

            self.patient_classes_ = classes
            self.patient_class_to_index_ = class_to_index

            if 1 in class_to_index:
                self.disease_class_index_ = class_to_index[1]
            elif "1" in class_to_index:
                self.disease_class_index_ = class_to_index["1"]
            else:
                self.disease_class_index_ = 1 if len(classes) > 1 else 0

            return "classification", {"y": y}, len(classes)

        else:
            raise ValueError("patient_task must be auto, cox, or classification.")

    def _prepare_cell_labels(self, cell_labels, cell_names, n_cells):
        if cell_labels is None:
            return None, None

        if isinstance(cell_labels, Mapping):
            if cell_names is None:
                raise ValueError("cell_names is required when cell_labels is dict.")

            cell_names = list(map(str, cell_names))
            raw = []

            for c in cell_names:
                if c in cell_labels:
                    raw.append(cell_labels[c])
                elif str(c) in cell_labels:
                    raw.append(cell_labels[str(c)])
                else:
                    raise KeyError(f"Cell {c} not found in cell_labels.")

            raw = np.asarray(raw)

        else:
            raw = np.asarray(cell_labels)

            if raw.ndim == 2:
                raw = np.argmax(raw, axis=1)

            if len(raw) != n_cells:
                raise ValueError("cell_labels length != number of cells.")

        classes = sorted(pd.unique(raw))
        class_to_index = {c: i for i, c in enumerate(classes)}
        y = np.asarray([class_to_index[v] for v in raw], dtype=np.int32)

        self.cell_classes_ = classes
        self.cell_class_to_index_ = class_to_index

        return y, len(classes)

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------

    def fit(
        self,
        sc_matrix,
        bulk_matrix,
        bulk_labels,
        sc_gene_names,
        bulk_gene_names,
        bulk_sample_names=None,
        cell_names=None,
        cell_labels=None,
    ):
        keras.backend.clear_session()

        sc_matrix = np.asarray(sc_matrix, dtype=np.float32)
        bulk_matrix = np.asarray(bulk_matrix, dtype=np.float32)

        if sc_matrix.ndim != 2 or bulk_matrix.ndim != 2:
            raise ValueError("sc_matrix and bulk_matrix must be 2D arrays.")

        n_cells = sc_matrix.shape[0]
        n_bulk = bulk_matrix.shape[0]

        if cell_names is None:
            cell_names = [f"cell_{i}" for i in range(n_cells)]
        else:
            cell_names = list(map(str, cell_names))

        if bulk_sample_names is not None:
            bulk_sample_names = list(map(str, bulk_sample_names))

        self.cell_names_ = cell_names
        self.bulk_sample_names_ = bulk_sample_names

        # Gene-name alignment
        sc_matrix, bulk_matrix, common_genes = align_genes(
            sc_matrix=sc_matrix,
            bulk_matrix=bulk_matrix,
            sc_gene_names=sc_gene_names,
            bulk_gene_names=bulk_gene_names,
        )

        self.common_genes_ = common_genes

        if self.log_transform:
            sc_matrix = np.log2(sc_matrix + 1.0)
            bulk_matrix = np.log2(bulk_matrix + 1.0)

        # Labels
        patient_task, bulk_label_data, n_patient_classes = self._prepare_bulk_labels(
            bulk_labels,
            bulk_sample_names,
            n_bulk,
        )
        self.patient_task_ = patient_task

        cell_y, n_cell_classes = self._prepare_cell_labels(
            cell_labels,
            cell_names,
            n_cells,
        )

        # Strict DEGAS-style feature selection
        self.selector_ = DEGASFeatureSelector(
            top_var_frac=self.top_var_frac,
            max_features=self.max_features,
            min_features=self.min_features,
            cox_l1_ratio=self.cox_l1_ratio,
            cox_n_alphas=self.cox_n_alphas,
            cox_alpha_min_ratio=self.cox_alpha_min_ratio,
            cox_max_iter=self.cox_max_iter,
            cox_tol=self.cox_tol,
            cox_top_n=self.cox_top_n,
        )

        selected_idx = self.selector_.select(
            sc_matrix=sc_matrix,
            bulk_matrix=bulk_matrix,
            patient_task=patient_task,
            bulk_label_data=bulk_label_data,
        )

        self.selected_gene_indices_ = selected_idx
        self.selected_genes_ = [self.common_genes_[i] for i in selected_idx]

        if self.verbose:
            ts_print(
                f"Common genes: {len(self.common_genes_)}",
                symbol="info",
            )
            ts_print(
                f"Selected genes: {len(self.selected_genes_)}",
                symbol="info",
            )

            if patient_task == "cox" and self.selector_.cox_nonzero_idx_ is not None:
                ts_print(
                    f"Elastic-net Cox nonzero genes: {len(self.selector_.cox_nonzero_idx_)}",
                    symbol="info",
                )

        sc_x = scale_expression(sc_matrix[:, selected_idx])
        bulk_x = scale_expression(bulk_matrix[:, selected_idx])

        sce = keras.losses.SparseCategoricalCrossentropy()

        self.models_ = []

        for b in range(self.n_bootstrap):
            seed = self.random_state + b
            np_rng = np.random.default_rng(seed)
            tf.random.set_seed(seed)

            model = DenseNetDEGAS(
                patient_task=patient_task,
                n_patient_classes=n_patient_classes,
                n_cell_classes=n_cell_classes,
                hidden_units=self.hidden_units,
                n_layers=self.n_layers,
                dropout_rate=self.dropout_rate,
                l2_base=self.l2_base,
            )

            optimizer = keras.optimizers.Adam(self.learning_rate)

            cell_pool = np_rng.choice(n_cells, size=n_cells, replace=True)
            bulk_pool = np_rng.choice(n_bulk, size=n_bulk, replace=True)

            for step in range(1, self.training_steps + 1):
                ci = np_rng.choice(
                    cell_pool,
                    size=min(self.cell_batch_size, len(cell_pool)),
                    replace=True,
                )
                bi = np_rng.choice(
                    bulk_pool,
                    size=min(self.bulk_batch_size, len(bulk_pool)),
                    replace=True,
                )

                x_cell = tf.convert_to_tensor(sc_x[ci], dtype=tf.float32)
                x_bulk = tf.convert_to_tensor(bulk_x[bi], dtype=tf.float32)

                with tf.GradientTape() as tape:
                    z_cell, _, cell_out = model(x_cell, training=True)
                    z_bulk, patient_out, _ = model(x_bulk, training=True)

                    if patient_task == "cox":
                        t = tf.convert_to_tensor(
                            bulk_label_data["time"][bi], dtype=tf.float32
                        )
                        s = tf.convert_to_tensor(
                            bulk_label_data["status"][bi], dtype=tf.float32
                        )
                        patient_loss = cox_ph_loss(t, s, patient_out)
                    else:
                        yb = tf.convert_to_tensor(
                            bulk_label_data["y"][bi], dtype=tf.int32
                        )
                        patient_loss = sce(yb, patient_out)

                    if cell_y is not None:
                        yc = tf.convert_to_tensor(cell_y[ci], dtype=tf.int32)
                        cell_loss = sce(yc, cell_out)
                    else:
                        cell_loss = tf.constant(0.0, dtype=tf.float32)

                    mmd_loss = mmd_rbf_loss(z_cell, z_bulk)

                    if model.losses:
                        reg_loss = tf.add_n(model.losses)
                    else:
                        reg_loss = tf.constant(0.0, dtype=tf.float32)

                    total_loss = (
                        self.lambda_patient * patient_loss
                        + self.lambda_cell * cell_loss
                        + self.lambda_mmd * mmd_loss
                        + self.lambda_l2 * reg_loss
                    )

                grads = tape.gradient(total_loss, model.trainable_variables)
                optimizer.apply_gradients(zip(grads, model.trainable_variables))

                if self.verbose and step % 500 == 0:
                    ts_print(
                        f"[bootstrap {b + 1}/{self.n_bootstrap}] "
                        f"step={step} "
                        f"total={float(total_loss.numpy()):.4f} "
                        f"patient={float(patient_loss.numpy()):.4f} "
                        f"cell={float(cell_loss.numpy()):.4f} "
                        f"mmd={float(mmd_loss.numpy()):.4f}",
                        symbol="debug",
                    )

            self.models_.append(model)

        return self

    # ------------------------------------------------------------------
    # Prediction
    # ------------------------------------------------------------------

    def _prepare_prediction_matrix(self, matrix, gene_names):
        if self.selected_genes_ is None:
            raise RuntimeError("Model has not been fitted.")

        matrix = np.asarray(matrix, dtype=np.float32)
        gene_names = np.asarray(gene_names).astype(str)

        gene_map = {}
        for i, g in enumerate(gene_names):
            if g not in gene_map:
                gene_map[g] = i

        missing = [g for g in self.selected_genes_ if g not in gene_map]

        if len(missing) > 0:
            raise ValueError(
                f"Prediction matrix is missing {len(missing)} selected genes. "
                f"Examples: {missing[:10]}"
            )

        idx = [gene_map[g] for g in self.selected_genes_]

        x = matrix[:, idx]

        if self.log_transform:
            x = np.log2(x + 1.0)

        return scale_expression(x)

    def predict_patient_head(self, matrix, gene_names, row_names=None, batch_size=2048):
        if not self.models_:
            raise RuntimeError("Model has not been fitted.")

        x = self._prepare_prediction_matrix(matrix, gene_names)
        n = x.shape[0]

        all_preds = []

        for model in self.models_:
            out = []

            for start in range(0, n, batch_size):
                end = min(start + batch_size, n)
                xb = tf.convert_to_tensor(x[start:end], dtype=tf.float32)
                _, patient_out, _ = model(xb, training=False)
                out.append(patient_out.numpy())

            all_preds.append(np.vstack(out))

        mean_pred = np.mean(all_preds, axis=0)

        if row_names is None:
            row_names = [f"sample_{i}" for i in range(n)]

        if self.patient_task_ == "cox":
            df = pd.DataFrame(
                {
                    "risk_score": mean_pred[:, 0],
                    "risk_association_centered": mean_pred[:, 0]
                    - np.mean(mean_pred[:, 0]),
                },
                index=row_names,
            )
        else:
            colnames = [f"class_{str(c)}_score" for c in self.patient_classes_]
            df = pd.DataFrame(mean_pred, columns=colnames, index=row_names)

            if mean_pred.shape[1] == 2:
                idx = self.disease_class_index_
                df["disease_score"] = mean_pred[:, idx]
                df["disease_association"] = 2.0 * mean_pred[:, idx] - 1.0

        return df

    def predict_cells(self, sc_matrix, sc_gene_names, cell_names=None, batch_size=2048):
        if cell_names is None:
            if (
                self.cell_names_ is not None
                and len(self.cell_names_) == np.asarray(sc_matrix).shape[0]
            ):
                cell_names = self.cell_names_
            else:
                cell_names = [
                    f"cell_{i}" for i in range(np.asarray(sc_matrix).shape[0])
                ]

        return self.predict_patient_head(
            matrix=sc_matrix,
            gene_names=sc_gene_names,
            row_names=list(map(str, cell_names)),
            batch_size=batch_size,
        )

    def predict_cell_head(self, matrix, gene_names, row_names=None, batch_size=2048):
        if not self.models_:
            raise RuntimeError("Model has not been fitted.")

        if self.cell_classes_ is None:
            raise RuntimeError("No cell label head was trained.")

        x = self._prepare_prediction_matrix(matrix, gene_names)
        n = x.shape[0]

        all_preds = []

        for model in self.models_:
            out = []

            for start in range(0, n, batch_size):
                end = min(start + batch_size, n)
                xb = tf.convert_to_tensor(x[start:end], dtype=tf.float32)
                _, _, cell_out = model(xb, training=False)
                out.append(cell_out.numpy())

            all_preds.append(np.vstack(out))

        mean_pred = np.mean(all_preds, axis=0)

        if row_names is None:
            row_names = [f"sample_{i}" for i in range(n)]

        colnames = [f"cell_class_{str(c)}_score" for c in self.cell_classes_]

        return pd.DataFrame(mean_pred, columns=colnames, index=row_names)

    def predict(
        self, matrix, gene_names, output="pat", row_names=None, batch_size=2048
    ):
        """
        Compatibility prediction.

        output:
            'pat' -> patient disease/survival head
            'sc'  -> cell-type head
        """

        if str(output).lower() in ["pat", "patient"]:
            return self.predict_patient_head(matrix, gene_names, row_names, batch_size)

        if str(output).lower() in ["sc", "cell"]:
            return self.predict_cell_head(matrix, gene_names, row_names, batch_size)

        raise ValueError("output must be 'pat' or 'sc'.")

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------

    def get_common_genes(self):
        return self.common_genes_

    def get_selected_genes(self):
        return self.selected_genes_

    def get_selected_gene_indices(self):
        return self.selected_gene_indices_

    def get_cox_nonzero_genes(self):
        if self.selector_ is None or self.selector_.cox_nonzero_idx_ is None:
            return None

        return [self.common_genes_[i] for i in self.selector_.cox_nonzero_idx_]

    def get_bulk_feature_scores(self):
        if self.selector_ is None:
            return None

        return pd.DataFrame(
            {
                "gene": self.common_genes_,
                "bulk_score": self.selector_.bulk_scores_,
                "sc_variance": self.selector_.sc_scores_,
            }
        )


# ----------------------------------------------------------------------
# Clean reticulate-facing API
# ----------------------------------------------------------------------


def train_degas(
    sc_matrix,
    bulk_matrix,
    bulk_labels,
    sc_gene_names,
    bulk_gene_names,
    bulk_sample_names=None,
    cell_names=None,
    cell_labels=None,
    patient_task="auto",
    **kwargs,
):
    model = DEGASTensorFlow(patient_task=patient_task, **kwargs)

    model.fit(
        sc_matrix=sc_matrix,
        bulk_matrix=bulk_matrix,
        bulk_labels=bulk_labels,
        sc_gene_names=sc_gene_names,
        bulk_gene_names=bulk_gene_names,
        bulk_sample_names=bulk_sample_names,
        cell_names=cell_names,
        cell_labels=cell_labels,
    )

    return model


def predict_degas(
    model,
    sc_matrix,
    sc_gene_names,
    cell_names=None,
    batch_size=2048,
):
    return model.predict_cells(
        sc_matrix=sc_matrix,
        sc_gene_names=sc_gene_names,
        cell_names=cell_names,
        batch_size=batch_size,
    )


def run_degas(
    sc_matrix,
    bulk_matrix,
    bulk_labels,
    sc_gene_names,
    bulk_gene_names,
    bulk_sample_names=None,
    cell_names=None,
    cell_labels=None,
    patient_task="auto",
    **kwargs,
):
    model = train_degas(
        sc_matrix=sc_matrix,
        bulk_matrix=bulk_matrix,
        bulk_labels=bulk_labels,
        sc_gene_names=sc_gene_names,
        bulk_gene_names=bulk_gene_names,
        bulk_sample_names=bulk_sample_names,
        cell_names=cell_names,
        cell_labels=cell_labels,
        patient_task=patient_task,
        **kwargs,
    )

    scores = model.predict_cells(
        sc_matrix=sc_matrix,
        sc_gene_names=sc_gene_names,
        cell_names=cell_names,
    )

    return scores


# ----------------------------------------------------------------------
# Original-name compatibility helpers
# ----------------------------------------------------------------------


def runCCMTLBag(
    scExp,
    scLab,
    patExp,
    patLab,
    tmpDir=None,
    model_type=None,
    architecture="DenseNet",
    FFdepth=3,
    Bagdepth=5,
    sc_gene_names=None,
    bulk_gene_names=None,
    bulk_sample_names=None,
    cell_names=None,
    patient_task="auto",
    **kwargs,
):
    """
    Compatibility-style wrapper.

    Unlike original DEGAS, this modern implementation does not write temporary
    TensorFlow graph files. tmpDir, model_type, and architecture are accepted
    only for interface familiarity.

    Required in this wrapper:
        sc_gene_names
        bulk_gene_names
    """

    if sc_gene_names is None or bulk_gene_names is None:
        raise ValueError(
            "Modern runCCMTLBag requires sc_gene_names and bulk_gene_names "
            "for strict gene alignment."
        )

    if model_type is not None:
        mt = str(model_type).lower()
        if "cox" in mt:
            patient_task = "cox"
        elif "class" in mt:
            patient_task = "classification"

    return train_degas(
        sc_matrix=scExp,
        bulk_matrix=patExp,
        bulk_labels=patLab,
        sc_gene_names=sc_gene_names,
        bulk_gene_names=bulk_gene_names,
        bulk_sample_names=bulk_sample_names,
        cell_names=cell_names,
        cell_labels=scLab,
        patient_task=patient_task,
        n_bootstrap=Bagdepth,
        n_layers=FFdepth,
        **kwargs,
    )


def predClassBag(
    ccModel,
    Exp,
    scORpat,
    gene_names,
    row_names=None,
    batch_size=2048,
):
    """
    Compatibility-style prediction wrapper.

    scORpat:
        'pat' -> disease/survival head
        'sc'  -> cell-type head
    """

    return ccModel.predict(
        matrix=Exp,
        gene_names=gene_names,
        output=scORpat,
        row_names=row_names,
        batch_size=batch_size,
    )
