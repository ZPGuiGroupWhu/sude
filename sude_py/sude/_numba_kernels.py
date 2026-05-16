import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist

try:
    from numba import get_num_threads, get_thread_id, njit, prange

    NUMBA_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised when numba is unavailable
    get_num_threads = None
    get_thread_id = None
    njit = None
    prange = None
    NUMBA_AVAILABLE = False


if NUMBA_AVAILABLE:

    @njit
    def _modified_distances_block_numba(
        distances,
        landmark_knn,
        rnn,
        agg_coef,
        start,
        stop,
        marker_size,
    ):
        n_samples = landmark_knn.shape[0]
        n_neighbors = landmark_knn.shape[1]
        block_size = stop - start
        output = np.empty((block_size, n_samples), np.float64)
        marker = np.zeros(marker_size, np.int64)
        snn = np.empty(n_samples, np.float64)
        tiny = np.finfo(np.float64).tiny

        for block_i in range(block_size):
            i = start + block_i
            for k in range(n_neighbors):
                marker[landmark_knn[i, k]] = 1

            max_snn = 0.0
            for j in range(n_samples):
                total = 0.0
                for k in range(n_neighbors):
                    neighbor = landmark_knn[j, k]
                    if marker[neighbor] == 1:
                        total += rnn[neighbor]
                snn[j] = total
                if total > max_snn:
                    max_snn = total

            if max_snn <= tiny:
                max_snn = tiny

            for j in range(n_samples):
                factor = 1.0 - snn[j] / max_snn
                output[block_i, j] = (factor**agg_coef) * distances[i, j]

            for k in range(n_neighbors):
                marker[landmark_knn[i, k]] = 0

        return output

    @njit
    def _build_neighbor_probability_rows_numba(
        distances,
        landmark_knn,
        rnn,
        k2,
        agg_coef,
        marker_size,
    ):
        n_samples = landmark_knn.shape[0]
        n_neighbors = landmark_knn.shape[1]
        rows = np.empty(n_samples * k2, np.int64)
        cols = np.empty(n_samples * k2, np.int64)
        values = np.empty(n_samples * k2, np.float64)
        marker = np.zeros(marker_size, np.int64)
        snn = np.empty(n_samples, np.float64)
        mod_dis = np.empty(n_samples, np.float64)
        selected = np.zeros(n_samples, np.int64)
        top_idx = np.empty(k2, np.int64)
        tiny = np.finfo(np.float64).tiny

        if k2 <= 0:
            return rows[:0], cols[:0], values[:0]

        out_pos = 0
        for i in range(n_samples):
            for k in range(n_neighbors):
                marker[landmark_knn[i, k]] = 1

            max_snn = 0.0
            for j in range(n_samples):
                total = 0.0
                for k in range(n_neighbors):
                    neighbor = landmark_knn[j, k]
                    if marker[neighbor] == 1:
                        total += rnn[neighbor]
                snn[j] = total
                if total > max_snn:
                    max_snn = total

            if max_snn <= tiny:
                max_snn = tiny

            for j in range(n_samples):
                factor = 1.0 - snn[j] / max_snn
                mod_dis[j] = (factor**agg_coef) * distances[i, j]

            mean_top = 0.0
            for k in range(k2):
                best_idx = -1
                best_value = np.inf
                for j in range(n_samples):
                    if selected[j] == 0:
                        value = mod_dis[j]
                        if value < best_value or (value == best_value and j < best_idx):
                            best_value = value
                            best_idx = j
                selected[best_idx] = 1
                top_idx[k] = best_idx
                mean_top += best_value
            mean_top /= k2
            denom = mean_top * mean_top
            if denom <= tiny:
                denom = tiny

            for k in range(k2):
                col = top_idx[k]
                value = np.exp(-0.5 * mod_dis[col] * mod_dis[col] / denom)
                rows[out_pos] = i
                cols[out_pos] = col
                values[out_pos] = value
                out_pos += 1

            for k in range(k2):
                selected[top_idx[k]] = 0

            for k in range(n_neighbors):
                marker[landmark_knn[i, k]] = 0

        return rows, cols, values

    @njit
    def _build_neighbor_probability_rows_streaming_numba(
        X_samp,
        landmark_knn,
        rnn,
        k2,
        agg_coef,
        marker_size,
    ):
        n_samples = landmark_knn.shape[0]
        n_features = X_samp.shape[1]
        n_neighbors = landmark_knn.shape[1]
        rows = np.empty(n_samples * k2, np.int64)
        cols = np.empty(n_samples * k2, np.int64)
        values = np.empty(n_samples * k2, np.float64)
        marker = np.zeros(marker_size, np.int64)
        snn = np.empty(n_samples, np.float64)
        mod_dis = np.empty(n_samples, np.float64)
        top_idx = np.empty(k2, np.int64)
        top_dist = np.empty(k2, np.float64)
        tiny = np.finfo(np.float64).tiny

        out_pos = 0
        for i in range(n_samples):
            for k in range(n_neighbors):
                marker[landmark_knn[i, k]] = 1

            max_snn = 0.0
            for j in range(n_samples):
                total = 0.0
                for k in range(n_neighbors):
                    neighbor = landmark_knn[j, k]
                    if marker[neighbor] == 1:
                        total += rnn[neighbor]
                snn[j] = total
                if total > max_snn:
                    max_snn = total

            if max_snn <= tiny:
                max_snn = tiny

            for k in range(k2):
                top_idx[k] = -1
                top_dist[k] = np.inf

            for j in range(n_samples):
                dist_sq = 0.0
                for feature in range(n_features):
                    delta = X_samp[i, feature] - X_samp[j, feature]
                    dist_sq += delta * delta
                factor = 1.0 - snn[j] / max_snn
                value = (factor**agg_coef) * np.sqrt(dist_sq)
                mod_dis[j] = value

                insert_at = -1
                for k in range(k2):
                    if value < top_dist[k] or (
                        value == top_dist[k] and (top_idx[k] == -1 or j < top_idx[k])
                    ):
                        insert_at = k
                        break
                if insert_at != -1:
                    for k in range(k2 - 1, insert_at, -1):
                        top_dist[k] = top_dist[k - 1]
                        top_idx[k] = top_idx[k - 1]
                    top_dist[insert_at] = value
                    top_idx[insert_at] = j

            mean_top = 0.0
            for k in range(k2):
                mean_top += top_dist[k]
            mean_top /= k2
            denom = mean_top * mean_top
            if denom <= tiny:
                denom = tiny

            for k in range(k2):
                col = top_idx[k]
                value = mod_dis[col]
                rows[out_pos] = i
                cols[out_pos] = col
                values[out_pos] = np.exp(-0.5 * value * value / denom)
                out_pos += 1

            for k in range(n_neighbors):
                marker[landmark_knn[i, k]] = 0

        return rows, cols, values

    @njit(parallel=True)
    def _fused_embedding_gradient_numba(indptr, indices, data, y):
        n_samples = y.shape[0]
        no_dims = y.shape[1]
        p_gradient = np.zeros((n_samples, no_dims), np.float64)
        q_gradient = np.zeros((n_samples, no_dims), np.float64)
        sum_q = 0.0

        for i in prange(n_samples):
            row_sum_q = 0.0
            for j in range(n_samples):
                dist_sq = 0.0
                for dim in range(no_dims):
                    delta = y[i, dim] - y[j, dim]
                    dist_sq += delta * delta
                q1 = 1.0 / (1.0 + np.log(1.0 + dist_sq))
                qq1 = 1.0 / (1.0 + dist_sq)
                coeff = 4.0 * q1 * q1 * qq1
                for dim in range(no_dims):
                    q_gradient[i, dim] += coeff * (y[i, dim] - y[j, dim])
                row_sum_q += q1
            sum_q += row_sum_q

            row_start = indptr[i]
            row_stop = indptr[i + 1]
            for ptr in range(row_start, row_stop):
                j = indices[ptr]
                p_value = data[ptr]
                dist_sq = 0.0
                for dim in range(no_dims):
                    delta = y[i, dim] - y[j, dim]
                    dist_sq += delta * delta
                q1 = 1.0 / (1.0 + np.log(1.0 + dist_sq))
                qq1 = 1.0 / (1.0 + dist_sq)
                coeff = 4.0 * p_value * q1 * qq1
                for dim in range(no_dims):
                    p_gradient[i, dim] += coeff * (y[i, dim] - y[j, dim])

        gradient = np.empty((n_samples, no_dims), np.float64)
        denom = sum_q - n_samples
        for i in prange(n_samples):
            for dim in range(no_dims):
                gradient[i, dim] = p_gradient[i, dim] - q_gradient[i, dim] / denom
        return gradient

    @njit(parallel=True)
    def _fused_embedding_gradient_2d_symmetric_numba(indptr, indices, data, y):
        n_samples = y.shape[0]
        p_gradient = np.zeros((n_samples, 2), np.float64)
        n_threads = get_num_threads()
        q_gradient_local = np.zeros((n_threads, n_samples, 2), np.float64)
        sum_q_local = np.zeros(n_threads, np.float64)

        for i in prange(n_samples):
            thread_id = get_thread_id()
            yi0 = y[i, 0]
            yi1 = y[i, 1]
            # Diagonal entries have zero gradient contribution but count in sum_q.
            sum_q_local[thread_id] += 1.0
            for j in range(i + 1, n_samples):
                dx = yi0 - y[j, 0]
                dy = yi1 - y[j, 1]
                dist_sq = dx * dx + dy * dy
                q1 = 1.0 / (1.0 + np.log(1.0 + dist_sq))
                qq1 = 1.0 / (1.0 + dist_sq)
                coeff = 4.0 * q1 * q1 * qq1
                gx = coeff * dx
                gy = coeff * dy
                q_gradient_local[thread_id, i, 0] += gx
                q_gradient_local[thread_id, i, 1] += gy
                q_gradient_local[thread_id, j, 0] -= gx
                q_gradient_local[thread_id, j, 1] -= gy
                sum_q_local[thread_id] += 2.0 * q1

            for ptr in range(indptr[i], indptr[i + 1]):
                j = indices[ptr]
                dx = yi0 - y[j, 0]
                dy = yi1 - y[j, 1]
                dist_sq = dx * dx + dy * dy
                q1 = 1.0 / (1.0 + np.log(1.0 + dist_sq))
                qq1 = 1.0 / (1.0 + dist_sq)
                coeff = 4.0 * data[ptr] * q1 * qq1
                p_gradient[i, 0] += coeff * dx
                p_gradient[i, 1] += coeff * dy

        q_gradient = np.zeros((n_samples, 2), np.float64)
        sum_q = 0.0
        for thread_id in range(n_threads):
            sum_q += sum_q_local[thread_id]
            for i in range(n_samples):
                q_gradient[i, 0] += q_gradient_local[thread_id, i, 0]
                q_gradient[i, 1] += q_gradient_local[thread_id, i, 1]

        gradient = np.empty((n_samples, 2), np.float64)
        denom = sum_q - n_samples
        for i in prange(n_samples):
            gradient[i, 0] = p_gradient[i, 0] - q_gradient[i, 0] / denom
            gradient[i, 1] = p_gradient[i, 1] - q_gradient[i, 1] / denom
        return gradient


def _build_neighbor_probability_rows_python(distances, landmark_knn, rnn, k2, agg_coef):
    n_samples = landmark_knn.shape[0]
    rows = []
    cols = []
    values = []
    knn_rnn_mat = rnn[landmark_knn]
    for i in range(n_samples):
        snn_id = np.isin(landmark_knn, landmark_knn[i]).astype(int)
        nn_id = np.where(np.max(snn_id, axis=1) == 1)[0]
        snn = np.zeros((1, n_samples))
        snn[:, nn_id] = np.sum(knn_rnn_mat[nn_id] * snn_id[nn_id], axis=1)
        mod_dis = (1 - snn / max(np.max(snn), np.finfo(float).tiny)) ** agg_coef * distances[
            i : i + 1,
            :,
        ]
        idx = np.argsort(mod_dis, axis=1)
        sort_dis = np.take_along_axis(mod_dis, idx, axis=1)
        mean_samp_dis_squared = np.square(np.mean(sort_dis[0, :k2]))
        values.extend(
            np.exp(
                -0.5
                * np.square(sort_dis[0, :k2])
                / np.maximum(mean_samp_dis_squared, np.finfo(float).tiny)
            )
        )
        rows.extend((i * np.ones((k2, 1))).flatten().tolist())
        cols.extend(idx[0, :k2])
    return np.asarray(rows), np.asarray(cols), np.asarray(values)


def _fused_embedding_gradient_python(probability, y):
    probability = probability.tocsr()
    y = np.asarray(y, dtype=np.float64)
    n_samples, no_dims = y.shape
    p_gradient = np.zeros((n_samples, no_dims))
    q_gradient = np.zeros((n_samples, no_dims))
    sum_q = 0.0

    for i in range(n_samples):
        for j in range(n_samples):
            diff = y[i] - y[j]
            dist_sq = float(diff @ diff)
            q1 = 1.0 / (1.0 + np.log(1.0 + dist_sq))
            qq1 = 1.0 / (1.0 + dist_sq)
            q_gradient[i] += 4.0 * q1 * q1 * qq1 * diff
            sum_q += q1

        for ptr in range(probability.indptr[i], probability.indptr[i + 1]):
            j = probability.indices[ptr]
            diff = y[i] - y[j]
            dist_sq = float(diff @ diff)
            q1 = 1.0 / (1.0 + np.log(1.0 + dist_sq))
            qq1 = 1.0 / (1.0 + dist_sq)
            p_gradient[i] += 4.0 * probability.data[ptr] * q1 * qq1 * diff

    return p_gradient - q_gradient / (sum_q - n_samples)


def fused_embedding_gradient(probability, y):
    probability = probability.tocsr()
    y = np.ascontiguousarray(y, dtype=np.float64)
    if NUMBA_AVAILABLE:
        return _fused_embedding_gradient_numba(
            probability.indptr,
            probability.indices,
            probability.data,
            y,
        )
    return _fused_embedding_gradient_python(probability, y)


def fused_embedding_gradient_2d_symmetric(probability, y):
    probability = probability.tocsr()
    y = np.ascontiguousarray(y, dtype=np.float64)
    if NUMBA_AVAILABLE:
        return _fused_embedding_gradient_2d_symmetric_numba(
            probability.indptr,
            probability.indices,
            probability.data,
            y,
        )
    return _fused_embedding_gradient_python(probability, y)


def build_neighbor_probability_rows(X_samp, landmark_knn, rnn, k2, agg_coef):
    distances = cdist(X_samp, X_samp)
    landmark_knn = np.asarray(landmark_knn, dtype=np.int64)
    rnn = np.asarray(rnn, dtype=np.float64)
    if NUMBA_AVAILABLE:
        marker_size = max(int(np.max(landmark_knn)) + 1, rnn.shape[0])
        n_samples = landmark_knn.shape[0]
        block_rows = min(256, n_samples)
        row_blocks = []
        col_blocks = []
        value_blocks = []
        for start in range(0, n_samples, block_rows):
            stop = min(start + block_rows, n_samples)
            mod_dis = _modified_distances_block_numba(
                distances,
                landmark_knn,
                rnn,
                float(agg_coef),
                start,
                stop,
                marker_size,
            )
            idx = np.argsort(mod_dis, axis=1)
            top_idx = idx[:, :k2]
            sort_dis = np.take_along_axis(mod_dis, idx, axis=1)[:, :k2]
            mean_squared = np.square(np.mean(sort_dis, axis=1))
            values = np.exp(
                -0.5
                * np.square(sort_dis)
                / np.maximum(mean_squared[:, np.newaxis], np.finfo(float).tiny)
            )
            row_blocks.append(np.repeat(np.arange(start, stop), k2))
            col_blocks.append(top_idx.reshape(-1))
            value_blocks.append(values.reshape(-1))
        return (
            np.concatenate(row_blocks),
            np.concatenate(col_blocks),
            np.concatenate(value_blocks),
        )
    return _build_neighbor_probability_rows_python(
        distances,
        landmark_knn,
        rnn,
        k2,
        agg_coef,
    )


def build_neighbor_probability_rows_streaming(X_samp, landmark_knn, rnn, k2, agg_coef):
    X_samp = np.ascontiguousarray(X_samp, dtype=np.float64)
    landmark_knn = np.asarray(landmark_knn, dtype=np.int64)
    rnn = np.asarray(rnn, dtype=np.float64)
    if NUMBA_AVAILABLE:
        marker_size = max(int(np.max(landmark_knn)) + 1, rnn.shape[0])
        return _build_neighbor_probability_rows_streaming_numba(
            X_samp,
            landmark_knn,
            rnn,
            int(k2),
            float(agg_coef),
            marker_size,
        )
    distances = cdist(X_samp, X_samp)
    return _build_neighbor_probability_rows_python(
        distances,
        landmark_knn,
        rnn,
        k2,
        agg_coef,
    )


def build_neighbor_probability_rows_sparse_snn(
    X_samp,
    landmark_knn,
    rnn,
    k2,
    agg_coef,
    block_rows=512,
):
    X_samp = np.asarray(X_samp, dtype=np.float64)
    landmark_knn = np.asarray(landmark_knn, dtype=np.int64)
    rnn = np.asarray(rnn, dtype=np.float64)
    n_samples = landmark_knn.shape[0]
    marker_size = max(int(np.max(landmark_knn)) + 1, rnn.shape[0])

    rows = np.repeat(np.arange(n_samples), landmark_knn.shape[1])
    cols = landmark_knn.reshape(-1)
    weights = rnn[cols]
    membership = csr_matrix(
        (weights, (rows, cols)),
        shape=(n_samples, marker_size),
    )
    snn = membership @ csr_matrix(
        (
            np.ones(landmark_knn.size),
            (cols, rows),
        ),
        shape=(marker_size, n_samples),
    )
    snn = snn.tocsr()

    row_blocks = []
    col_blocks = []
    value_blocks = []
    for start in range(0, n_samples, block_rows):
        stop = min(start + block_rows, n_samples)
        distances = cdist(X_samp[start:stop], X_samp)
        snn_block = snn[start:stop].toarray()
        max_snn = np.maximum(np.max(snn_block, axis=1), np.finfo(float).tiny)
        modified = ((1.0 - snn_block / max_snn[:, np.newaxis]) ** agg_coef) * distances
        if k2 < n_samples:
            top_idx_unsorted = np.argpartition(modified, kth=k2 - 1, axis=1)[:, :k2]
            top_dis_unsorted = np.take_along_axis(modified, top_idx_unsorted, axis=1)
            order = np.argsort(top_dis_unsorted, axis=1)
            top_idx = np.take_along_axis(top_idx_unsorted, order, axis=1)
            top_dis = np.take_along_axis(top_dis_unsorted, order, axis=1)
        else:
            top_idx = np.argsort(modified, axis=1)
            top_dis = np.take_along_axis(modified, top_idx, axis=1)
        mean_squared = np.square(np.mean(top_dis, axis=1))
        values = np.exp(
            -0.5
            * np.square(top_dis)
            / np.maximum(mean_squared[:, np.newaxis], np.finfo(float).tiny)
        )
        row_blocks.append(np.repeat(np.arange(start, stop), k2))
        col_blocks.append(top_idx.reshape(-1))
        value_blocks.append(values.reshape(-1))

    return (
        np.concatenate(row_blocks),
        np.concatenate(col_blocks),
        np.concatenate(value_blocks),
    )
