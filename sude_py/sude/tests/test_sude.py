import unittest

import numpy as np
from scipy.spatial.distance import pdist
from sklearn.datasets import load_iris
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from sude import SUDE, sude


class TestSUDE(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.X = load_iris(return_X_y=True)[0][:20]
        cls.X_query = np.vstack([cls.X[:4], cls.X[:2]])

    def assert_embeddings_equivalent(self, embedding_a, embedding_b):
        self.assertEqual(embedding_a.shape, embedding_b.shape)
        np.testing.assert_allclose(
            pdist(embedding_a),
            pdist(embedding_b),
            atol=1e-10,
            rtol=1e-10,
        )

    def test_function_matches_estimator_fit_transform(self):
        model = SUDE(
            n_components=2,
            n_neighbors=0,
            init="pca",
            max_iter=2,
        )
        embedding_from_estimator = model.fit_transform(self.X)
        embedding_from_function = sude(
            self.X,
            no_dims=2,
            k1=0,
            initialize="pca",
            T_epoch=2,
        )

        np.testing.assert_allclose(
            embedding_from_estimator,
            embedding_from_function,
        )

    def test_fit_sets_expected_attributes(self):
        model = SUDE(
            n_components=2,
            n_neighbors=10,
            init="pca",
            max_iter=1,
        ).fit(self.X)

        self.assertEqual(model.embedding_.shape, (20, 2))
        self.assertEqual(model.Y_landmarks_.shape[1], 2)
        self.assertEqual(model.n_features_in_, self.X.shape[1])
        self.assertGreater(model.n_landmarks_, 0)
        self.assertEqual(model.get_feature_names_out().tolist(), ["sude0", "sude1"])

    def test_fit_caches_landmark_neighbor_index(self):
        model = SUDE(
            n_components=2,
            n_neighbors=10,
            init="pca",
            max_iter=1,
        ).fit(self.X)

        self.assertTrue(hasattr(model, "landmark_nn_"))

    def test_returns_expected_shape(self):
        embedding = sude(self.X, no_dims=2, k1=0, initialize="pca", T_epoch=2)

        self.assertEqual(embedding.shape, (20, 2))
        self.assertTrue(np.isfinite(embedding).all())

    def test_duplicate_rows_keep_duplicate_embeddings(self):
        X = np.vstack([self.X[:10], self.X[:3]])

        embedding = sude(X, no_dims=2, k1=0, initialize="pca", T_epoch=1)

        np.testing.assert_allclose(embedding[0], embedding[10])
        np.testing.assert_allclose(embedding[1], embedding[11])
        np.testing.assert_allclose(embedding[2], embedding[12])

    def test_estimator_supported_initializers(self):
        X = self.X[:15]

        for init in ("spectral", "le", "pca", "mds"):
            with self.subTest(init=init):
                embedding = SUDE(
                    n_components=2,
                    n_neighbors=0,
                    init=init,
                    max_iter=1,
                ).fit_transform(X)
                self.assertEqual(embedding.shape, (15, 2))
                self.assertTrue(np.isfinite(embedding).all())

    def test_transform_embeds_new_samples(self):
        model = SUDE(
            n_components=2,
            n_neighbors=10,
            init="pca",
            max_iter=1,
        ).fit(self.X)

        embedding = model.transform(self.X_query)

        self.assertEqual(embedding.shape, (6, 2))
        self.assertTrue(np.isfinite(embedding).all())
        np.testing.assert_allclose(embedding[0], embedding[4])
        np.testing.assert_allclose(embedding[1], embedding[5])

    def test_pipeline_compatibility(self):
        pipeline = make_pipeline(
            StandardScaler(),
            SUDE(
                n_components=2,
                n_neighbors=0,
                init="pca",
                max_iter=1,
                normalize=False,
            ),
        )
        embedding = pipeline.fit_transform(self.X)

        self.assertEqual(embedding.shape, (20, 2))
        self.assertTrue(np.isfinite(embedding).all())

    def test_function_keeps_legacy_initializer_names(self):
        embedding = sude(
            self.X[:15],
            no_dims=2,
            k1=0,
            initialize="le",
            T_epoch=1,
        )
        estimator_embedding = SUDE(
            n_components=2,
            n_neighbors=0,
            init="spectral",
            max_iter=1,
        ).fit_transform(self.X[:15])

        self.assertEqual(embedding.shape, (15, 2))
        self.assert_embeddings_equivalent(embedding, estimator_embedding)

    def test_function_supported_initializers(self):
        X = self.X[:15]

        for initialize in ("le", "pca", "mds"):
            with self.subTest(initialize=initialize):
                embedding = sude(
                    X,
                    no_dims=2,
                    k1=0,
                    initialize=initialize,
                    T_epoch=1,
                )
                self.assertEqual(embedding.shape, (15, 2))
                self.assertTrue(np.isfinite(embedding).all())

    def test_large_mode_runs(self):
        embedding = sude(
            self.X,
            no_dims=2,
            k1=0,
            initialize="pca",
            large=True,
            T_epoch=1,
        )

        self.assertEqual(embedding.shape, (20, 2))
        self.assertTrue(np.isfinite(embedding).all())


if __name__ == "__main__":
    unittest.main()
