import unittest

import numpy as np
from sklearn.datasets import load_iris

from sude import SUDE, sude
from sude._sude import _compute_landmarks
from sude.clle import clle, clle_batch
from sude.pps import pps


def pps_reference(knn, rnn, order):
    id_samp = []
    id_sort = sorted(range(len(rnn)), key=lambda k: rnn[k], reverse=True)
    while len(id_sort) > 0:
        id_samp.append(id_sort[0])
        rm_pts = [id_sort[0]]
        for _ in range(order):
            rm_pts.extend(knn[rm_pts].flatten().tolist())
        rm_pts = set(rm_pts)
        rm_id = np.where(np.isin(id_sort, list(rm_pts)))[0]
        id_sort = [id_sort[i] for i in range(len(id_sort)) if i not in rm_id]
    return id_samp


def compute_landmarks_reference(X, n_components, n_neighbors):
    from sklearn.neighbors import NearestNeighbors

    from sude.init_pca import init_pca

    n_samples, n_features = X.shape
    if n_neighbors == 0:
        return [], [], list(range(n_samples))

    if n_samples > 5000 and n_features > 50:
        reduced = init_pca(X, n_components, 0.8)
        get_knn = (
            NearestNeighbors(n_neighbors=n_neighbors + 1)
            .fit(reduced)
            .kneighbors(reduced, return_distance=False)
        )
    else:
        get_knn = (
            NearestNeighbors(n_neighbors=n_neighbors + 1)
            .fit(X)
            .kneighbors(X, return_distance=False)
        )
    _, rnn = np.unique(get_knn, return_counts=True)
    id_samp = pps_reference(get_knn, rnn, 1)
    return get_knn, rnn, id_samp


def clle_reference(X_samp, Y_samp, X_i, N_dis):
    n = X_samp.shape[0]
    S = (X_samp - X_i) @ (X_samp - X_i).transpose()
    if np.abs(np.linalg.det(S)) <= np.finfo(float).eps:
        S = S + (0.1**2 / n) * np.trace(S) * np.eye(n)
    W = (np.linalg.inv(S) @ np.ones((n, 1))) / (
        np.ones((1, n)) @ np.linalg.inv(S) @ np.ones((n, 1))
    )
    Y_0 = W.transpose() @ Y_samp
    dd = np.sqrt((Y_samp[0] - Y_0) @ (Y_samp[0] - Y_0).transpose())
    if dd != 0:
        Y_i = Y_samp[0] + N_dis * (Y_0 - Y_samp[0]) / dd
    else:
        Y_i = Y_samp[0]

    return Y_i


class TestReferenceEquivalence(unittest.TestCase):
    def test_pps_matches_reference(self):
        knn = np.array(
            [
                [0, 1, 2],
                [1, 0, 2],
                [2, 1, 3],
                [3, 2, 4],
                [4, 3, 2],
            ]
        )
        rnn = np.array([3, 2, 4, 1, 1])

        self.assertEqual(pps(knn, rnn, 1), pps_reference(knn, rnn, 1))

    def test_compute_landmarks_matches_reference(self):
        X = load_iris(return_X_y=True)[0][:40]

        actual_knn, actual_rnn, actual_ids = _compute_landmarks(X, 2, 5)
        expected_knn, expected_rnn, expected_ids = compute_landmarks_reference(X, 2, 5)

        np.testing.assert_array_equal(actual_knn, expected_knn)
        np.testing.assert_array_equal(actual_rnn, expected_rnn)
        self.assertEqual(actual_ids, expected_ids)

    def test_clle_matches_reference(self):
        X_samp = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        Y_samp = np.array([[0.0, 0.0], [0.5, 0.0], [0.0, 0.5]])
        X_i = np.array([0.25, 0.25])

        actual = clle(X_samp, Y_samp, X_i, 0.5)
        expected = clle_reference(X_samp, Y_samp, X_i, 0.5)

        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)

    def test_clle_batch_matches_scalar_reference(self):
        rng = np.random.RandomState(11)
        n_queries = 17
        n_neighbors = 4
        n_features = 6
        n_components = 3
        X_neighbors = rng.rand(n_queries, n_neighbors, n_features)
        Y_neighbors = rng.rand(n_queries, n_neighbors, n_components)
        X_query = rng.rand(n_queries, n_features)
        n_dis = rng.rand(n_queries)

        actual = clle_batch(X_neighbors, Y_neighbors, X_query, n_dis)
        expected = np.vstack(
            [
                clle_reference(
                    X_neighbors[i],
                    Y_neighbors[i],
                    X_query[i],
                    n_dis[i],
                )
                for i in range(n_queries)
            ]
        )

        np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10)

    def test_sude_function_and_estimator_match_for_pca(self):
        rng = np.random.RandomState(1)
        X = rng.rand(40, 6)

        estimator = SUDE(
            n_components=2,
            n_neighbors=10,
            init="pca",
            max_iter=2,
        ).fit_transform(X)
        legacy = sude(X, no_dims=2, k1=10, initialize="pca", T_epoch=2)

        np.testing.assert_allclose(estimator, legacy, rtol=1e-12, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
