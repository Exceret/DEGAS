import os

os.environ.setdefault("KERAS_BACKEND", "tensorflow")
import numpy as np
import pandas as pd
import tensorflow as tf
import keras

# -----------------------------
# Loss functions
# -----------------------------


def cox_ph_loss(time, status, pred_risk):
    """
    Negative Cox partial log-likelihood.
    time: shape [batch]
    status: 1 event, 0 censored
    pred_risk: model risk score, shape [batch] or [batch, 1]
    """
    time = tf.cast(tf.reshape(time, [-1]), tf.float32)
    status = tf.cast(tf.reshape(status, [-1]), tf.float32)
    pred_risk = tf.cast(tf.reshape(pred_risk, [-1]), tf.float32)

    order = tf.argsort(time, direction="DESCENDING")
    risk_sorted = tf.gather(pred_risk, order)
    status_sorted = tf.gather(status, order)

    log_cumsum_risk = tf.math.cumulative_logsumexp(risk_sorted)
    loss = -tf.reduce_sum((risk_sorted - log_cumsum_risk) * status_sorted)

    n_events = tf.reduce_sum(status_sorted) + 1e-8
    return loss / n_events


def _pairwise_sq_dists(x, y):
    x_norm = tf.reduce_sum(tf.square(x), axis=1, keepdims=True)
    y_norm = tf.reduce_sum(tf.square(y), axis=1, keepdims=True)
    d = x_norm - 2.0 * tf.matmul(x, y, transpose_b=True) + tf.transpose(y_norm)
    return tf.maximum(d, 0.0)


def mmd_rbf_loss(x, y, sigmas=(1.0, 2.0, 4.0, 8.0, 16.0)):
    """
    Maximum Mean Discrepancy with multi-scale RBF kernels.
    """
    x = tf.cast(x, tf.float32)
    y = tf.cast(y, tf.float32)

    d_xx = _pairwise_sq_dists(x, x)
    d_yy = _pairwise_sq_dists(y, y)
    d_xy = _pairwise_sq_dists(x, y)

    loss = 0.0
    for sigma in sigmas:
        beta = 1.0 / (2.0 * sigma**2)
        k_xx = tf.exp(-beta * d_xx)
        k_yy = tf.exp(-beta * d_yy)
        k_xy = tf.exp(-beta * d_xy)

        loss += tf.reduce_mean(k_xx) + tf.reduce_mean(k_yy) - 2.0 * tf.reduce_mean(k_xy)

    return loss / float(len(sigmas))


# -----------------------------
# DEGAS model
# -----------------------------


class _DenseNetDEGAS(keras.Model):
    def __init__(
        self,
        patient_task,
        n_patient_classes=2,
        n_cell_classes=None,
        hidden_units=50,
        n_layers=3,
        dropout_rate=0.5,
        l2_base=1e-4,
    ):
        super().__init__()

        self.patient_task = patient_task
        self.n_cell_classes = n_cell_classes

        reg = keras.regularizers.l2(l2_base) if l2_base > 0 else None

        self.hidden_layers_ = [
            keras.layers.Dense(
                hidden_units,
                activation="sigmoid",
                kernel_regularizer=reg,
                bias_regularizer=reg,
            )
            for _ in range(n_layers)
        ]

        self.dropout = keras.layers.Dropout(dropout_rate)

        if patient_task == "cox":
            self.patient_head = keras.layers.Dense(
                1,
                activation="sigmoid",
                kernel_regularizer=reg,
                bias_regularizer=reg,
            )
        elif patient_task == "classification":
            self.patient_head = keras.layers.Dense(
                n_patient_classes,
                activation="softmax",
                kernel_regularizer=reg,
                bias_regularizer=reg,
            )
        else:
            raise ValueError("patient_task must be 'cox' or 'classification'.")

        if n_cell_classes is not None:
            self.cell_head = keras.layers.Dense(
                n_cell_classes,
                activation="softmax",
                kernel_regularizer=reg,
                bias_regularizer=reg,
            )
        else:
            self.cell_head = None

    def encode(self, x, training=False):
        """
        DenseNet-style encoder.
        Each hidden layer receives concatenated previous representations.
        Final hidden layer is used as latent space for MMD.
        """
        concat_h = x
        last_h = None

        for layer in self.hidden_layers_:
            h = layer(concat_h)
            h = self.dropout(h, training=training)
            concat_h = tf.concat([concat_h, h], axis=1)
            last_h = h

        return last_h

    def call(self, x, training=False):
        z = self.encode(x, training=training)
        patient_out = self.patient_head(z)

        if self.cell_head is not None:
            cell_out = self.cell_head(z)
        else:
            cell_out = None

        return z, patient_out, cell_out


# -----------------------------
# Main API
# -----------------------------


class DEGASTensorFlow:
    """
    TensorFlow implementation of DEGAS-like transfer learning.

    Main use:
        model = DEGASTensorFlow(patient_task="cox")
        model.fit(sc_matrix, bulk_matrix, bulk_labels)
        cell_scores = model.predict_cells(sc_matrix)
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
    ):
        self.patient_task = patient_task
        self.n_bootstrap = n_bootstrap
        self.training_steps = training_steps
        self.cell_batch_size = cell_batch_size
        self.bulk_batch_size = bulk_batch_size
        self.hidden_units = hidden_units
        self.n_layers = n_layers
        self.dropout_rate = dropout_rate
        self.lambda_cell = lambda_cell
        self.lambda_patient = lambda_patient
        self.lambda_mmd = lambda_mmd
        self.lambda_l2 = lambda_l2
        self.l2_base = l2_base
        self.learning_rate = learning_rate
        self.top_var_frac = top_var_frac
        self.max_features = max_features
        self.min_features = min_features
        self.random_state = random_state
        self.verbose = verbose

        self.models_ = []
        self.feature_idx_ = None
        self.patient_task_ = None
        self.cell_names_ = None
        self.bulk_sample_names_ = None
        self.cell_label_mapping_ = None

    # ---------- preprocessing ----------

    def _select_features(self, sc_matrix, bulk_matrix):
        n_genes = sc_matrix.shape[1]

        sc_var = np.nanvar(sc_matrix, axis=0)
        bulk_var = np.nanvar(bulk_matrix, axis=0)

        k = int(max(self.min_features, n_genes * self.top_var_frac))
        k = min(k, n_genes)

        sc_top = np.argpartition(-sc_var, k - 1)[:k]
        bulk_top = np.argpartition(-bulk_var, k - 1)[:k]

        inter = np.intersect1d(sc_top, bulk_top)

        if len(inter) < self.min_features:
            combined_var = sc_var + bulk_var
            k2 = min(max(self.min_features, len(inter)), n_genes, self.max_features)
            inter = np.argpartition(-combined_var, k2 - 1)[:k2]

        if self.max_features is not None and len(inter) > self.max_features:
            combined_var = sc_var[inter] + bulk_var[inter]
            order = np.argsort(-combined_var)[: self.max_features]
            inter = inter[order]

        return np.sort(inter)

    @staticmethod
    def _scale_expression(x):
        """
        Paper-style:
        1. sample-wise z-score
        2. sample-wise min-max scaling to [0,1]
        """
        x = np.asarray(x, dtype=np.float32)
        x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)

        mean = x.mean(axis=1, keepdims=True)
        std = x.std(axis=1, keepdims=True)
        std[std == 0] = 1.0

        z = (x - mean) / std

        z_min = z.min(axis=1, keepdims=True)
        z_max = z.max(axis=1, keepdims=True)
        denom = z_max - z_min
        denom[denom == 0] = 1.0

        return ((z - z_min) / denom).astype(np.float32)

    # ---------- label preparation ----------

    def _infer_patient_task(self, bulk_labels):
        if self.patient_task != "auto":
            return self.patient_task

        if isinstance(bulk_labels, pd.DataFrame):
            return "cox"
        elif isinstance(bulk_labels, dict):
            return "classification"
        else:
            raise ValueError(
                "Cannot infer patient_task. Use patient_task='cox' or 'classification'."
            )

    def _prepare_bulk_labels(self, bulk_labels, bulk_sample_names, n_bulk):
        task = self._infer_patient_task(bulk_labels)

        if task == "cox":
            if not isinstance(bulk_labels, pd.DataFrame):
                raise ValueError("For Cox task, bulk_labels must be pd.DataFrame.")

            df = bulk_labels.copy()
            colmap = {c.lower(): c for c in df.columns}

            if "time" not in colmap or "status" not in colmap:
                raise ValueError(
                    "Cox label DataFrame must contain columns: time, status."
                )

            if bulk_sample_names is not None:
                bulk_sample_names = list(bulk_sample_names)
                if set(bulk_sample_names).issubset(set(df.index)):
                    df = df.loc[bulk_sample_names]
                elif len(df) != n_bulk:
                    raise ValueError(
                        "bulk_sample_names not found in Cox label index, "
                        "and label length != bulk samples."
                    )

            if len(df) != n_bulk:
                raise ValueError("Cox label rows must equal number of bulk samples.")

            time = df[colmap["time"]].values.astype(np.float32)
            status = df[colmap["status"]].values.astype(np.float32)

            return task, {"time": time, "status": status}, 1

        elif task == "classification":
            if isinstance(bulk_labels, dict):
                if bulk_sample_names is None:
                    raise ValueError(
                        "For dict labels, bulk_sample_names must be provided "
                        "and match rows of bulk_matrix."
                    )

                y = []
                for name in bulk_sample_names:
                    if name in bulk_labels:
                        y.append(bulk_labels[name])
                    elif str(name) in bulk_labels:
                        y.append(bulk_labels[str(name)])
                    else:
                        raise KeyError(f"Sample {name} not found in bulk_labels dict.")
                y = np.asarray(y, dtype=np.int32)

            else:
                y = np.asarray(bulk_labels, dtype=np.int32)
                if len(y) != n_bulk:
                    raise ValueError("Classification labels length != bulk samples.")

            if not np.all(np.isin(y, [0, 1])):
                raise ValueError("This implementation expects binary labels 0/1.")

            return task, {"y": y}, 2

        else:
            raise ValueError("patient_task must be 'cox' or 'classification'.")

    def _prepare_cell_labels(self, cell_labels, cell_names, n_cells):
        """
        Optional. If not provided, DEGAS becomes BlankCox or BlankClass.
        """
        if cell_labels is None:
            return None, None

        if isinstance(cell_labels, dict):
            if cell_names is None:
                raise ValueError("cell_names required when cell_labels is dict.")

            raw = []
            for name in cell_names:
                if name in cell_labels:
                    raw.append(cell_labels[name])
                elif str(name) in cell_labels:
                    raw.append(cell_labels[str(name)])
                else:
                    raise KeyError(f"Cell {name} not found in cell_labels dict.")
            raw = np.asarray(raw)
        else:
            raw = np.asarray(cell_labels)
            if len(raw) != n_cells:
                raise ValueError("cell_labels length != number of cells.")

        classes = sorted(pd.unique(raw))
        mapping = {c: i for i, c in enumerate(classes)}
        y = np.asarray([mapping[v] for v in raw], dtype=np.int32)

        self.cell_label_mapping_ = mapping

        return y, len(classes)

    # ---------- fit ----------

    def fit(
        self,
        sc_matrix,
        bulk_matrix,
        bulk_labels,
        bulk_sample_names=None,
        cell_names=None,
        cell_labels=None,
    ):
        """
        Parameters
        ----------
        sc_matrix : np.ndarray
            cells × genes
        bulk_matrix : np.ndarray
            samples × genes
        bulk_labels :
            Cox: pd.DataFrame with columns ['time', 'status']
            Classification: dict {sample_name: 0/1}
        bulk_sample_names : list-like
            Required when bulk_labels is dict.
        cell_names : list-like
            Optional cell names.
        cell_labels :
            Optional cell subtype labels. If None, no cell classification head is trained.
        """

        keras.backend.clear_session()

        sc_matrix = np.asarray(sc_matrix, dtype=np.float32)
        bulk_matrix = np.asarray(bulk_matrix, dtype=np.float32)

        if sc_matrix.ndim != 2 or bulk_matrix.ndim != 2:
            raise ValueError("sc_matrix and bulk_matrix must be 2D numpy arrays.")

        if sc_matrix.shape[1] != bulk_matrix.shape[1]:
            raise ValueError("sc_matrix and bulk_matrix must have same gene columns.")

        n_cells, n_genes = sc_matrix.shape
        n_bulk = bulk_matrix.shape[0]

        if cell_names is None:
            cell_names = [f"cell_{i}" for i in range(n_cells)]
        if bulk_sample_names is not None:
            bulk_sample_names = list(bulk_sample_names)

        self.cell_names_ = list(cell_names)
        self.bulk_sample_names_ = bulk_sample_names

        patient_task, prepared_bulk_labels, n_patient_classes = (
            self._prepare_bulk_labels(bulk_labels, bulk_sample_names, n_bulk)
        )
        self.patient_task_ = patient_task

        prepared_cell_labels, n_cell_classes = self._prepare_cell_labels(
            cell_labels, cell_names, n_cells
        )

        # feature selection
        self.feature_idx_ = self._select_features(sc_matrix, bulk_matrix)

        sc_x = self._scale_expression(sc_matrix[:, self.feature_idx_])
        bulk_x = self._scale_expression(bulk_matrix[:, self.feature_idx_])

        # n_features = sc_x.shape[1]

        # rng_global = np.random.default_rng(self.random_state)
        self.models_ = []

        sce = keras.losses.SparseCategoricalCrossentropy()

        for b in range(self.n_bootstrap):
            seed = self.random_state + b
            tf.random.set_seed(seed)
            rng = np.random.default_rng(seed)

            model = _DenseNetDEGAS(
                patient_task=patient_task,
                n_patient_classes=n_patient_classes,
                n_cell_classes=n_cell_classes,
                hidden_units=self.hidden_units,
                n_layers=self.n_layers,
                dropout_rate=self.dropout_rate,
                l2_base=self.l2_base,
            )

            optimizer = keras.optimizers.Adam(self.learning_rate)

            # bootstrap pools
            cell_pool = rng.choice(n_cells, size=n_cells, replace=True)
            bulk_pool = rng.choice(n_bulk, size=n_bulk, replace=True)

            for step in range(1, self.training_steps + 1):
                ci = rng.choice(
                    cell_pool,
                    size=min(self.cell_batch_size, len(cell_pool)),
                    replace=True,
                )
                bi = rng.choice(
                    bulk_pool,
                    size=min(self.bulk_batch_size, len(bulk_pool)),
                    replace=True,
                )

                x_cell = tf.convert_to_tensor(sc_x[ci], dtype=tf.float32)
                x_bulk = tf.convert_to_tensor(bulk_x[bi], dtype=tf.float32)

                with tf.GradientTape() as tape:
                    z_cell, _, cell_out = model(x_cell, training=True)
                    z_bulk, patient_out, _ = model(x_bulk, training=True)

                    # patient loss
                    if patient_task == "cox":
                        t = tf.convert_to_tensor(
                            prepared_bulk_labels["time"][bi], dtype=tf.float32
                        )
                        s = tf.convert_to_tensor(
                            prepared_bulk_labels["status"][bi], dtype=tf.float32
                        )
                        patient_loss = cox_ph_loss(t, s, patient_out)
                    else:
                        yb = tf.convert_to_tensor(
                            prepared_bulk_labels["y"][bi], dtype=tf.int32
                        )
                        patient_loss = sce(yb, patient_out)

                    # optional cell classification loss
                    if prepared_cell_labels is not None:
                        yc = tf.convert_to_tensor(
                            prepared_cell_labels[ci], dtype=tf.int32
                        )
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
                    print(
                        f"[bootstrap {b+1}/{self.n_bootstrap}] "
                        f"step={step} "
                        f"total={total_loss.numpy():.4f} "
                        f"patient={patient_loss.numpy():.4f} "
                        f"cell={cell_loss.numpy():.4f} "
                        f"mmd={mmd_loss.numpy():.4f}"
                    )

            self.models_.append(model)

        return self

    # ---------- predict ----------

    def predict_cells(self, sc_matrix, cell_names=None, batch_size=2048):
        """
        Predict patient disease attribute scores for each single cell.

        Returns
        -------
        pd.DataFrame
        """
        if not self.models_:
            raise RuntimeError("Model has not been fitted.")

        sc_matrix = np.asarray(sc_matrix, dtype=np.float32)
        x = self._scale_expression(sc_matrix[:, self.feature_idx_])

        n_cells = x.shape[0]

        if cell_names is None:
            if self.cell_names_ is not None and len(self.cell_names_) == n_cells:
                cell_names = self.cell_names_
            else:
                cell_names = [f"cell_{i}" for i in range(n_cells)]

        all_preds = []

        for model in self.models_:
            preds = []

            for start in range(0, n_cells, batch_size):
                end = min(start + batch_size, n_cells)
                xb = tf.convert_to_tensor(x[start:end], dtype=tf.float32)
                _, patient_out, _ = model(xb, training=False)
                preds.append(patient_out.numpy())

            preds = np.vstack(preds)
            all_preds.append(preds)

        mean_pred = np.mean(all_preds, axis=0)

        if self.patient_task_ == "cox":
            df = pd.DataFrame(
                {
                    "risk_score": mean_pred[:, 0],
                    "risk_association_centered": mean_pred[:, 0]
                    - np.mean(mean_pred[:, 0]),
                },
                index=cell_names,
            )

        else:
            df = pd.DataFrame(
                mean_pred,
                columns=[f"class_{i}_score" for i in range(mean_pred.shape[1])],
                index=cell_names,
            )

            if mean_pred.shape[1] == 2:
                df["disease_score"] = mean_pred[:, 1]
                df["disease_association"] = 2.0 * mean_pred[:, 1] - 1.0

        return df

    def get_selected_gene_indices(self):
        return self.feature_idx_


def main() -> None:
    # ----------------------------------------------------------------
    # 1. Load data
    # ----------------------------------------------------------------
    data_dir = "/data/home/yyx/Project/SigBridgeR_methods/DEGAS/DEGAS_res"

    # Single-cell expression: cells × genes
    sc_mat = pd.read_csv(os.path.join(data_dir, "scExp.csv"), index_col=0)
    print(f"[INFO] scExp loaded: {sc_mat.shape[0]} cells × {sc_mat.shape[1]} genes")

    # Bulk/patient expression: samples × genes (must have same gene columns)
    bulk_mat = pd.read_csv(os.path.join(data_dir, "patExp.csv"), index_col=0)
    print(
        f"[INFO] patExp loaded: {bulk_mat.shape[0]} samples × {bulk_mat.shape[1]} genes"
    )

    # Patient phenotype (Cox: time, status)
    phenotype = pd.read_csv(os.path.join(data_dir, "phenotype.csv"), index_col=False)
    print(
        f"[INFO] phenotype loaded: {phenotype.shape[0]} samples, columns={list(phenotype.columns)}"
    )

    # Ensure sc_mat and bulk_mat share the same gene set and order
    common_genes = sc_mat.columns.intersection(bulk_mat.columns).sort_values()
    print(f"[INFO] Common genes: {len(common_genes)}")
    sc_mat = sc_mat[common_genes]
    bulk_mat = bulk_mat[common_genes]

    # Convert to numpy
    sc_matrix = sc_mat.values.astype(np.float32)
    bulk_matrix = bulk_mat.values.astype(np.float32)

    # ----------------------------------------------------------------
    # 2. Initialize and fit DEGAS model
    # ----------------------------------------------------------------
    # phenotype has 'time' and 'status' columns -> Cox task
    model = DEGASTensorFlow(
        patient_task="cox",  # Cox regression for survival data
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
    )

    print("\n[INFO] Starting DEGAS training...")
    model.fit(
        sc_matrix=sc_matrix,
        bulk_matrix=bulk_matrix,
        bulk_labels=phenotype,  # pd.DataFrame with 'time' and 'status'
        bulk_sample_names=bulk_mat.index.tolist(),
    )

    # ----------------------------------------------------------------
    # 3. Predict cell-level risk scores
    # ----------------------------------------------------------------
    print("\n[INFO] Predicting cell-level risk scores...")
    cell_scores = model.predict_cells(sc_matrix)

    # ----------------------------------------------------------------
    # 4. Save results
    # ----------------------------------------------------------------
    output_dir = os.path.join(data_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    # Save cell scores
    cell_scores.to_csv(os.path.join(output_dir, "cell_risk_scores.csv"))
    print(f"[INFO] Cell risk scores saved to {output_dir}/cell_risk_scores.csv")

    # Save selected gene indices
    selected_genes = model.get_selected_gene_indices()
    gene_names = common_genes[selected_genes]
    pd.DataFrame({"gene": gene_names, "index": selected_genes}).to_csv(
        os.path.join(output_dir, "selected_genes.csv"), index=False
    )
    print(
        f"[INFO] Selected genes ({len(selected_genes)}) saved to {output_dir}/selected_genes.csv"
    )

    # Print summary
    print("\n" + "=" * 60)
    print("DEGAS training completed!")
    print(f"  Patient task: {model.patient_task_}")
    print(f"  Number of bootstrap models: {len(model.models_)}")
    print(f"  Selected features: {len(selected_genes)}")
    print(f"  Cell scores shape: {cell_scores.shape}")
    print(f"  Cell scores columns: {list(cell_scores.columns)}")
    print(f"  Top 5 risk scores:\n{cell_scores.head()}")
    print("=" * 60)


if __name__ == "__main__":
    main()
