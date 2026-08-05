import tensorflow as tf


def cox_ph_loss(time, status, pred_risk):
    """
    Negative Cox partial log-likelihood.

    Parameters
    ----------
    time : tensor, shape [n]
        Follow-up or survival time.
    status : tensor, shape [n]
        Event indicator. 1 = event, 0 = censored.
    pred_risk : tensor, shape [n] or [n, 1]
        Neural-network Cox risk score.

    Returns
    -------
    Tensor scalar.
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
    Maximum Mean Discrepancy using multi-scale RBF kernels.
    Computed between latent cell representation and latent bulk representation.
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
