import os

os.environ.setdefault("KERAS_BACKEND", "tensorflow")

import tensorflow as tf
import keras


class DenseNetDEGAS(keras.Model):
    """
    DEGAS DenseNet-style multitask model.

    Shared encoder:
        cells and bulk samples use same input genes and same hidden representation.

    Heads:
        patient_head: Cox or classification.
        cell_head: optional cell-type classification.
    """

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

        reg = keras.regularizers.l2(l2_base) if l2_base and l2_base > 0 else None

        self.hidden_layers_ = [
            keras.layers.Dense(
                hidden_units,
                activation="sigmoid",
                kernel_regularizer=reg,
                bias_regularizer=reg,
                name=f"degas_hidden_{i + 1}",
            )
            for i in range(n_layers)
        ]

        self.dropout = keras.layers.Dropout(dropout_rate)

        if patient_task == "cox":
            self.patient_head = keras.layers.Dense(
                1,
                activation="sigmoid",
                kernel_regularizer=reg,
                bias_regularizer=reg,
                name="patient_cox_head",
            )
        elif patient_task == "classification":
            self.patient_head = keras.layers.Dense(
                n_patient_classes,
                activation="softmax",
                kernel_regularizer=reg,
                bias_regularizer=reg,
                name="patient_class_head",
            )
        else:
            raise ValueError("patient_task must be 'cox' or 'classification'.")

        if n_cell_classes is not None:
            self.cell_head = keras.layers.Dense(
                n_cell_classes,
                activation="softmax",
                kernel_regularizer=reg,
                bias_regularizer=reg,
                name="cell_class_head",
            )
        else:
            self.cell_head = None

    def encode(self, x, training=False):
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
