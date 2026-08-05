import numpy as np
from scipy import stats
from sklearn.preprocessing import StandardScaler

try:
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sksurv.util import Surv
except Exception:
    CoxnetSurvivalAnalysis = None
    Surv = None

from .preprocessing import rank_normalize


class DEGASFeatureSelector:
    """
    DEGAS-style feature selection.

    scRNA:
        high-variance genes.

    bulk classification:
        binary: Welch t-test.
        multiclass: one-way ANOVA.

    bulk Cox:
        strict elastic-net Cox using scikit-survival CoxnetSurvivalAnalysis.

    Final:
        intersection of sc-selected and bulk-selected genes.
    """

    def __init__(
        self,
        top_var_frac=0.2,
        max_features=1000,
        min_features=50,
        cox_l1_ratio=0.5,
        cox_n_alphas=100,
        cox_alpha_min_ratio=0.01,
        cox_max_iter=100000,
        cox_tol=1e-7,
        cox_top_n=None,
    ):
        self.top_var_frac = top_var_frac
        self.max_features = max_features
        self.min_features = min_features

        self.cox_l1_ratio = cox_l1_ratio
        self.cox_n_alphas = cox_n_alphas
        self.cox_alpha_min_ratio = cox_alpha_min_ratio
        self.cox_max_iter = cox_max_iter
        self.cox_tol = cox_tol
        self.cox_top_n = cox_top_n

        self.sc_scores_ = None
        self.bulk_scores_ = None
        self.sc_selected_ = None
        self.bulk_selected_ = None
        self.selected_idx_ = None
        self.cox_nonzero_idx_ = None

    def select(self, sc_matrix, bulk_matrix, patient_task, bulk_label_data):
        sc_matrix = np.asarray(sc_matrix, dtype=np.float32)
        bulk_matrix = np.asarray(bulk_matrix, dtype=np.float32)

        if sc_matrix.shape[1] != bulk_matrix.shape[1]:
            raise ValueError("sc_matrix and bulk_matrix must have same gene count.")

        n_genes = sc_matrix.shape[1]

        # scRNA high-variance genes
        sc_scores = np.nanvar(sc_matrix, axis=0)
        k_sc = int(max(self.min_features, n_genes * self.top_var_frac))
        k_sc = min(k_sc, n_genes)
        sc_selected = np.argpartition(-sc_scores, k_sc - 1)[:k_sc]

        # bulk disease-attribute-associated genes
        k_bulk = int(max(self.min_features, n_genes * self.top_var_frac))
        k_bulk = min(k_bulk, n_genes)

        if patient_task == "classification":
            bulk_scores = self._classification_scores(
                bulk_matrix,
                bulk_label_data["y"],
            )
            bulk_selected = np.argpartition(-bulk_scores, k_bulk - 1)[:k_bulk]

        elif patient_task == "cox":
            bulk_selected, bulk_scores = self._select_bulk_genes_by_elastic_net_cox(
                bulk_matrix=bulk_matrix,
                time=bulk_label_data["time"],
                status=bulk_label_data["status"],
                n_select=k_bulk,
            )

        else:
            raise ValueError("patient_task must be 'classification' or 'cox'.")

        inter = np.intersect1d(sc_selected, bulk_selected)

        # If intersection is too small, fill by combined rank.
        # Cox information still comes strictly from elastic-net Cox coefficients.
        if len(inter) < self.min_features:
            combined = rank_normalize(sc_scores) + rank_normalize(bulk_scores)

            k_final = min(max(self.min_features, len(inter)), n_genes)
            if self.max_features is not None:
                k_final = min(k_final, self.max_features)

            inter = np.argpartition(-combined, k_final - 1)[:k_final]

        # If too many, truncate by combined score
        if self.max_features is not None and len(inter) > self.max_features:
            combined = rank_normalize(sc_scores) + rank_normalize(bulk_scores)
            local_score = combined[inter]
            order = np.argsort(-local_score)[: self.max_features]
            inter = inter[order]

        self.sc_scores_ = sc_scores
        self.bulk_scores_ = bulk_scores
        self.sc_selected_ = np.sort(sc_selected)
        self.bulk_selected_ = np.sort(bulk_selected)
        self.selected_idx_ = np.sort(inter)

        return self.selected_idx_

    @staticmethod
    def _classification_scores(x, y):
        x = np.asarray(x, dtype=np.float32)
        y = np.asarray(y)

        n_genes = x.shape[1]
        classes = np.unique(y)
        scores = np.zeros(n_genes, dtype=np.float32)

        if len(classes) < 2:
            raise ValueError("Classification requires at least two classes.")

        if len(classes) == 2:
            g0 = x[y == classes[0]]
            g1 = x[y == classes[1]]

            for j in range(n_genes):
                try:
                    stat, _ = stats.ttest_ind(
                        g0[:, j],
                        g1[:, j],
                        equal_var=False,
                        nan_policy="omit",
                    )
                    scores[j] = 0.0 if np.isnan(stat) else abs(stat)
                except Exception:
                    scores[j] = 0.0

        else:
            groups = [x[y == c] for c in classes]

            for j in range(n_genes):
                try:
                    vals = [g[:, j] for g in groups if g.shape[0] > 1]
                    if len(vals) < 2:
                        scores[j] = 0.0
                    else:
                        stat, _ = stats.f_oneway(*vals)
                        scores[j] = 0.0 if np.isnan(stat) else stat
                except Exception:
                    scores[j] = 0.0

        return np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)

    def _elastic_net_cox_scores(self, x, time, status):
        """
        Strict elastic-net Cox feature selection.

        Importance score:
            max(abs(coef)) across Coxnet regularization path.
        """

        if CoxnetSurvivalAnalysis is None or Surv is None:
            raise ImportError(
                "scikit-survival is required for strict elastic-net Cox feature selection. "
                "Install with: pip install scikit-survival "
                "or conda install -c conda-forge scikit-survival"
            )

        x = np.asarray(x, dtype=np.float32)
        time = np.asarray(time, dtype=np.float64)
        status = np.asarray(status).astype(bool)

        if len(time) != x.shape[0] or len(status) != x.shape[0]:
            raise ValueError("time/status length must equal number of bulk samples.")

        if np.sum(status) < 2:
            raise ValueError("Elastic-net Cox requires at least two observed events.")

        gene_std = np.nanstd(x, axis=0)
        valid = gene_std > 1e-8

        scores = np.zeros(x.shape[1], dtype=np.float32)

        if np.sum(valid) == 0:
            raise ValueError("All genes have near-zero variance in bulk matrix.")

        x_valid = x[:, valid]
        x_valid = np.nan_to_num(x_valid, nan=0.0, posinf=0.0, neginf=0.0)

        scaler = StandardScaler()
        x_scaled = scaler.fit_transform(x_valid)

        y_surv = Surv.from_arrays(
            event=status,
            time=time,
        )

        model = CoxnetSurvivalAnalysis(
            l1_ratio=self.cox_l1_ratio,
            n_alphas=self.cox_n_alphas,
            alpha_min_ratio=self.cox_alpha_min_ratio,
            max_iter=self.cox_max_iter,
            tol=self.cox_tol,
            fit_baseline_model=False,
        )

        model.fit(x_scaled, y_surv)

        coef_path = model.coef_

        if coef_path.shape[0] == x_scaled.shape[1]:
            coef_abs = np.abs(coef_path)
        elif coef_path.shape[1] == x_scaled.shape[1]:
            coef_abs = np.abs(coef_path.T)
        else:
            raise RuntimeError(
                f"Unexpected Coxnet coef_ shape: {coef_path.shape}, "
                f"n_features={x_scaled.shape[1]}"
            )

        valid_scores = np.max(coef_abs, axis=1)

        if np.all(valid_scores == 0):
            raise RuntimeError(
                "Elastic-net Cox selected zero genes. "
                "Try lowering cox_alpha_min_ratio, lowering cox_l1_ratio, "
                "or increasing sample size/events."
            )

        scores[valid] = valid_scores.astype(np.float32)
        return np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)

    def _select_bulk_genes_by_elastic_net_cox(
        self, bulk_matrix, time, status, n_select
    ):
        scores = self._elastic_net_cox_scores(
            bulk_matrix,
            time,
            status,
        )

        if self.cox_top_n is not None:
            n_select = min(n_select, int(self.cox_top_n))

        nonzero = np.where(scores > 0)[0]
        self.cox_nonzero_idx_ = nonzero

        if len(nonzero) >= n_select:
            local_scores = scores[nonzero]
            order = np.argsort(-local_scores)[:n_select]
            selected = nonzero[order]
        else:
            selected = np.argsort(-scores)[:n_select]

        return np.sort(selected), scores
