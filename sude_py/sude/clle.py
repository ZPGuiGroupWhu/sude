import numpy as np


def _regularize_singular_gram(S):
    n = S.shape[-1]
    if S.ndim == 2:
        if np.abs(np.linalg.det(S)) <= np.finfo(float).eps:
            S = S + (0.1**2 / n) * np.trace(S) * np.eye(n)
        return S

    singular = np.abs(np.linalg.det(S)) <= np.finfo(float).eps
    if np.any(singular):
        eye = np.eye(n, dtype=S.dtype)
        trace = np.trace(S[singular], axis1=1, axis2=2)
        S = np.array(S, copy=True)
        S[singular] += ((0.1**2 / n) * trace)[:, np.newaxis, np.newaxis] * eye
    return S


def clle(X_samp, Y_samp, X_i, N_dis):
    """
    Constrained Locally Linear Embedding (CLLE)
    This function returns representation of point X_i.

    Parameters are:

    'X_samp'      - High-dimensional features of KNN of point X_i. Each row denotes an observation.
    'Y_samp'      - Low-dimensional embeddings of KNN of point X_i.
    'X_i'         - Current non-landmark point.
    'N_dis'       - Distance between point X_i and its nearest neighbor in lower-dimensional space.

    """
    n = X_samp.shape[0]
    S = (X_samp - X_i) @ (X_samp - X_i).transpose()
    S = _regularize_singular_gram(S)
    ones = np.ones((n, 1))
    W = np.linalg.solve(S, ones)
    W = W / (ones.transpose() @ W)
    Y_0 = W.transpose() @ Y_samp
    dd = np.sqrt((Y_samp[0] - Y_0) @ (Y_samp[0] - Y_0).transpose())
    if dd != 0:
        Y_i = Y_samp[0] + N_dis * (Y_0 - Y_samp[0]) / dd
    else:
        Y_i = Y_samp[0]

    return Y_i


def clle_batch(X_samp, Y_samp, X_i, N_dis):
    """
    Batched CLLE for multiple non-landmark points.

    Each row in X_samp/Y_samp contains the nearest landmark neighbors for the
    corresponding row in X_i. The math matches clle(), but solves all tiny
    linear systems in one NumPy call and avoids explicit matrix inversion.
    """
    X_samp = np.asarray(X_samp, dtype=np.float64)
    Y_samp = np.asarray(Y_samp, dtype=np.float64)
    X_i = np.asarray(X_i, dtype=np.float64)
    N_dis = np.asarray(N_dis, dtype=np.float64)

    if X_i.shape[0] == 0:
        return np.empty((0, Y_samp.shape[2]), dtype=np.float64)

    n = X_samp.shape[1]
    diff = X_samp - X_i[:, np.newaxis, :]
    S = diff @ np.swapaxes(diff, 1, 2)
    S = _regularize_singular_gram(S)

    ones = np.ones((X_i.shape[0], n, 1), dtype=np.float64)
    W = np.linalg.solve(S, ones)
    W = W / np.sum(W, axis=1, keepdims=True)
    Y_0 = np.sum(W * Y_samp, axis=1)

    nearest_y = Y_samp[:, 0, :]
    delta = Y_0 - nearest_y
    dd = np.linalg.norm(delta, axis=1)
    Y_i = nearest_y.copy()
    nonzero = dd != 0
    Y_i[nonzero] = (
        nearest_y[nonzero]
        + N_dis[nonzero, np.newaxis] * delta[nonzero] / dd[nonzero, np.newaxis]
    )
    return Y_i
