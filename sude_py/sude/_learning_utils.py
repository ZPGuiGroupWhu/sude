import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh


def sorted_distances_and_indices(distances):
    idx = np.argsort(distances, axis=1)
    return np.take_along_axis(distances, idx, axis=1), idx


def apply_graph_gradient(pro_mat_y, y):
    column_sums = np.sum(pro_mat_y, axis=0)
    return column_sums[:, np.newaxis] * y - pro_mat_y @ y


def init_le_from_probability(probability, no_dims):
    degree = np.asarray(probability.sum(axis=0)).ravel()
    degree_sqrt = np.sqrt(degree)
    laplacian = diags(degree) - probability
    laplacian = laplacian.multiply(degree_sqrt[:, np.newaxis])
    laplacian = laplacian.multiply(degree_sqrt[np.newaxis, :])
    eigenvalues, eigenvectors = eigsh(laplacian, k=no_dims + 1, which="SM")
    smallest_indices = np.argsort(np.abs(eigenvalues))
    return eigenvectors[:, smallest_indices[1:]]
