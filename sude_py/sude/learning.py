from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist
from ._learning_utils import init_le_from_probability, sorted_distances_and_indices
from ._numba_kernels import (
    NUMBA_AVAILABLE,
    build_neighbor_probability_rows_sparse_snn,
    build_neighbor_probability_rows,
    fused_embedding_gradient,
    fused_embedding_gradient_2d_symmetric,
)
from .init_pca import init_pca
from .mds import mds
from .pca import pca
import math
import numpy as np
import time


def _adaptive_k2(n_samples):
    if n_samples < 9:
        return n_samples
    if n_samples > 1000:
        return int(np.ceil(np.log2(n_samples)) + 18)
    if n_samples > 50:
        return int(np.ceil(0.02 * n_samples)) + 8
    return 9


def _auto_block_size(n_samples, memory_budget_mb):
    if memory_budget_mb is None:
        return n_samples
    arrays_needed = 5
    bytes_per_float = np.dtype(np.float64).itemsize
    block_size = int(memory_budget_mb * 1024**2 / (n_samples * bytes_per_float * arrays_needed))
    return max(1, min(n_samples, block_size))


def _build_probability_matrix(
    X_samp,
    k1,
    get_knn,
    rnn,
    id_samp,
    no_dims,
    agg_coef,
    k2,
    memory_budget_mb=None,
):
    n_samples, n_features = X_samp.shape
    if k1 > 0:
        if NUMBA_AVAILABLE:
            landmark_knn = get_knn[id_samp]
            if memory_budget_mb is None:
                row, col, values = build_neighbor_probability_rows(
                    X_samp,
                    landmark_knn,
                    rnn,
                    k2,
                    agg_coef,
                )
            else:
                row, col, values = build_neighbor_probability_rows_sparse_snn(
                    X_samp,
                    landmark_knn,
                    rnn,
                    k2,
                    agg_coef,
                    block_rows=_auto_block_size(n_samples, memory_budget_mb),
                )
            return csr_matrix((values, (row, col)), shape=(n_samples, n_samples))

        row = []
        col = []
        values = []
        landmark_knn = get_knn[id_samp]
        knn_rnn_mat = rnn[landmark_knn]
        for i in range(n_samples):
            snn_id = np.isin(landmark_knn, landmark_knn[i]).astype(int)
            nn_id = np.where(np.max(snn_id, axis=1) == 1)[0]
            snn = np.zeros((1, n_samples))
            snn[:, nn_id] = np.sum(knn_rnn_mat[nn_id] * snn_id[nn_id], axis=1)
            mod_dis = (
                (1 - snn / max(np.max(snn), np.finfo(float).tiny)) ** agg_coef
                * cdist(X_samp[i : i + 1, :], X_samp)
            )
            sort_dis, idx = sorted_distances_and_indices(mod_dis)
            mean_samp_dis_squared = np.square(np.mean(sort_dis[0, :k2]))
            values.extend(
                np.exp(
                    -0.5
                    * np.square(sort_dis[0, :k2])
                    / np.maximum(mean_samp_dis_squared, np.finfo(float).tiny)
                )
            )
            row.extend((i * np.ones((k2, 1))).flatten().tolist())
            col.extend(idx[0, :k2])
        return csr_matrix((values, (row, col)), shape=(n_samples, n_samples))

    if n_samples > 5000 and n_features > 50:
        reduced = init_pca(X_samp, no_dims, 0.8)
        samp_dis, samp_knn = NearestNeighbors(n_neighbors=k2).fit(reduced).kneighbors(reduced)
    else:
        samp_dis, samp_knn = NearestNeighbors(n_neighbors=k2).fit(X_samp).kneighbors(X_samp)
    mean_samp_dis_squared = np.square(np.mean(samp_dis, axis=1))
    values = np.exp(
        -0.5
        * np.square(samp_dis)
        / np.maximum(mean_samp_dis_squared[:, np.newaxis], np.finfo(float).tiny)
    )
    return csr_matrix(
        (
            values.flatten(),
            ([i for i in range(n_samples) for _ in range(k2)], samp_knn.flatten()),
        ),
        shape=(n_samples, n_samples),
    )


def _initialize_embedding(X_samp, probability, no_dims, initialize):
    if initialize == "le":
        return init_le_from_probability(probability, no_dims)
    if initialize == "pca":
        return pca(X_samp, no_dims)
    if initialize == "mds":
        return mds(X_samp, no_dims)
    raise ValueError("initialize must be one of {'le', 'pca', 'mds'}")


def _dense_gradient(probability_dense, y, alpha, momentum, previous_gradient):
    distances = cdist(y, y) ** 2
    q1 = 1 / (1 + np.log(1 + distances))
    qq1 = 1 / (1 + distances)
    q = q1 / (np.sum(q1) - y.shape[0])
    pro_mat_y = 4 * (probability_dense - q) * q1 * qq1
    gradient = np.sum(pro_mat_y, axis=0)[:, np.newaxis] * y - pro_mat_y @ y
    return y - alpha * (gradient + momentum * previous_gradient), gradient


def _block_gradient(probability, y, alpha, momentum, previous_gradient, block_size):
    n_samples, no_dims = y.shape
    p_gradient = np.zeros((n_samples, no_dims))
    q_gradient = np.zeros((n_samples, no_dims))
    sum_q = 0.0
    for start in range(0, n_samples, block_size):
        stop = min(start + block_size, n_samples)
        idx = np.arange(start, stop)
        distances = cdist(y[idx], y) ** 2
        q1 = 1 / (1 + np.log(1 + distances))
        qq1 = 1 / (1 + distances)
        p_mat = -4 * probability[idx, :].multiply(q1).multiply(qq1).toarray()
        q_mat = -4 * q1**2 * qq1
        len_blk = idx.shape[0]
        id_pq = np.column_stack((np.arange(len_blk), start + np.arange(len_blk)))
        p_mat[id_pq[:, 0], id_pq[:, 1]] = p_mat[id_pq[:, 0], id_pq[:, 1]] - np.sum(
            p_mat,
            axis=1,
        )
        q_mat[id_pq[:, 0], id_pq[:, 1]] = q_mat[id_pq[:, 0], id_pq[:, 1]] - np.sum(
            q_mat,
            axis=1,
        )
        p_gradient[idx] = p_mat @ y
        q_gradient[idx] = q_mat @ y
        sum_q += np.sum(q1)
    gradient = p_gradient - q_gradient / (sum_q - n_samples)
    return y - alpha * (gradient + momentum * previous_gradient), gradient


def _fused_gradient(probability, y, alpha, momentum, previous_gradient):
    if y.shape[1] == 2:
        gradient = fused_embedding_gradient_2d_symmetric(probability, y)
    else:
        gradient = fused_embedding_gradient(probability, y)
    return y - alpha * (gradient + momentum * previous_gradient), gradient


def learning(
    X_samp,
    k1,
    get_knn,
    rnn,
    id_samp,
    no_dims,
    initialize,
    agg_coef,
    T_epoch,
    memory_budget_mb=None,
    return_profile=False,
):
    profile = {}
    total_start = time.perf_counter()
    n_samples, _ = X_samp.shape
    k2 = _adaptive_k2(n_samples)
    probability_start = time.perf_counter()
    probability = _build_probability_matrix(
        X_samp,
        k1,
        get_knn,
        rnn,
        id_samp,
        no_dims,
        agg_coef,
        k2,
        memory_budget_mb=memory_budget_mb,
    )
    probability = (probability + probability.transpose()) / 2
    profile["probability_seconds"] = time.perf_counter() - probability_start
    initialize_start = time.perf_counter()
    y = _initialize_embedding(X_samp, probability, no_dims, initialize)
    profile["initialize_seconds"] = time.perf_counter() - initialize_start
    normalize_start = time.perf_counter()
    probability = probability / (probability.sum() - n_samples)
    profile["normalize_probability_seconds"] = time.perf_counter() - normalize_start

    block_size = _auto_block_size(n_samples, memory_budget_mb)
    use_fused_gradient = NUMBA_AVAILABLE
    use_dense_gradient = (not use_fused_gradient) and block_size >= n_samples
    probability_dense = probability.toarray() if use_dense_gradient else None

    max_alpha = 2.5 * n_samples
    min_alpha = 2 * n_samples
    warm_step = 10
    previous_gradient = np.zeros((n_samples, no_dims))
    epoch = 1
    gradient_seconds = 0.0
    while epoch <= T_epoch:
        if epoch <= warm_step:
            alpha = max_alpha
        else:
            alpha = min_alpha + 0.5 * (max_alpha - min_alpha) * (
                1 + np.cos(np.pi * ((epoch - warm_step) / (T_epoch - warm_step)))
            )
        momentum = (epoch - 1) / (epoch + 2)
        gradient_start = time.perf_counter()
        if use_fused_gradient:
            y, gradient = _fused_gradient(probability, y, alpha, momentum, previous_gradient)
        elif use_dense_gradient:
            y, gradient = _dense_gradient(probability_dense, y, alpha, momentum, previous_gradient)
        else:
            y, gradient = _block_gradient(
                probability,
                y,
                alpha,
                momentum,
                previous_gradient,
                block_size,
            )
        gradient_seconds += time.perf_counter() - gradient_start
        previous_gradient = gradient
        epoch += 1

    profile["gradient_seconds"] = gradient_seconds
    profile["total_seconds"] = time.perf_counter() - total_start
    profile["n_samples"] = n_samples
    profile["k2"] = k2
    profile["used_fused_gradient"] = bool(use_fused_gradient)
    print(str(epoch - 1) + " epochs have been computed!")
    if return_profile:
        return y, k2, profile
    return y, k2


def memory_budget_for_large(n_samples):
    no_blocks = math.ceil(n_samples / 3000)
    block_size = math.ceil(n_samples / no_blocks)
    arrays_needed = 5
    bytes_per_float = np.dtype(np.float64).itemsize
    return block_size * n_samples * bytes_per_float * arrays_needed / 1024**2
