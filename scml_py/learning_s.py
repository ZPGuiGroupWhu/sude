from sklearn.neighbors import NearestNeighbors
from scipy.sparse import diags
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
from scipy.spatial.distance import cdist
from init_pca import init_pca
from pca import pca
from mds import mds
import numpy as np


def learning_s(X_samp, k1, get_knn, rnn, id_samp, no_dims, initialize, agg_coef, T_epoch, T_vcc):
    """
    This function returns representation of the landmarks in the lower-dimensional space and the number of nearest
    neighbors of landmarks. It computes the gradient using the entire probability matrix P and Q.

    """
    # Obtain size and dimension of landmarks
    N, dim = X_samp.shape

    # Compute the number of nearest neighbors of landmarks adaptively
    if N < 9:
        k2 = N
    else:
        if N > 1000:
            k2 = int(np.ceil(np.log2(N)) + 18)
        elif N > 50:
            k2 = int(np.ceil(0.02 * N)) + 8
        else:
            k2 = 9

    # Compute high-dimensional probability matrix P
    if k1 > 0:
        # Compute SNN matrix of landmarks
        SNN = np.zeros((N, N))
        knn_rnn_mat = rnn[get_knn[id_samp]]
        for i in range(N):
            snn_id = np.isin(get_knn[id_samp], get_knn[id_samp[i]]).astype(int)
            nn_id = np.where(np.max(snn_id, axis=1) == 1)[0]
            SNN[i, nn_id] = np.sum(knn_rnn_mat[nn_id] * snn_id[nn_id], axis=1)
            SNN[i] = SNN[i] / max(np.max(SNN[i]), np.finfo(float).tiny)
        # Compute the modified distance matrix
        Dis = (1 - SNN) ** agg_coef * cdist(X_samp, X_samp)
        P = np.zeros((N, N))
        sort_dis = np.sort(Dis, axis=1)
        idx = np.argsort(Dis, axis=1)
        for i in range(N):
            P[i, idx[i, :k2]] = np.exp(
                -0.5 * np.square(sort_dis[i, :k2]) / np.maximum(np.square(np.mean(sort_dis[i, :k2])),
                                                                np.finfo(float).tiny))
    else:
        if N > 5000 and dim > 50:
            xx = init_pca(X_samp, no_dims, 0.8)
            samp_dis, samp_knn = NearestNeighbors(n_neighbors=k2).fit(xx).kneighbors(xx)
        else:
            samp_dis, samp_knn = NearestNeighbors(n_neighbors=k2).fit(X_samp).kneighbors(X_samp)
        mean_samp_dis_squared = np.square(np.mean(samp_dis, axis=1))
        Pval = np.exp(-0.5 * np.square(samp_dis) / np.maximum(mean_samp_dis_squared[:, np.newaxis], np.finfo(float).tiny))
        P = csr_matrix((Pval.flatten(), ([i for i in range(N) for _ in range(k2)], samp_knn.flatten())), shape=(N, N)).toarray()
    # Symmetrize matrix P
    P = (P + P.transpose()) / 2

    # Initialize embedding Y of landmarks
    if initialize == 'le':
        Dg = diags(np.array(np.sum(P, axis=0)))
        L = np.sqrt(Dg) @ (Dg - P) @ np.sqrt(Dg)
        eigenvalues, eigenvectors = eigs(L, k=no_dims + 1, which='SM')
        smallest_indices = np.argsort(np.abs(eigenvalues))
        Y = np.real(eigenvectors[:, smallest_indices[1:]])
        del Dg, L
    elif initialize == 'pca':
        Y = pca(X_samp, no_dims)
    elif initialize == 'mds':
        Y = mds(X_samp, no_dims)

    # Normalize matrix P
    P = P / (np.sum(P) - N)

    # Compute matrix Q
    D = cdist(Y, Y) ** 2
    Q1 = 1 / (1 + np.log(1 + D))
    QQ1 = 1 / (1 + D)
    del D
    Q = Q1 / (np.sum(Q1) - N)

    # Initialization
    alpha = 2.5 * N
    preGrad = np.zeros((N, no_dims))
    KL = -P * np.log(Q)
    KL[np.isnan(KL)] = 0
    cost = np.sum(KL)
    vcc = 1
    epoch = 1
    length = 3
    while epoch <= T_epoch and vcc > T_vcc:
        # Compute gradient
        ProMatY = 4 * (P - Q) * Q1 * QQ1
        grad = (np.diag(np.sum(ProMatY, axis=0)) - ProMatY) @ Y
        # Update embedding Y
        Y = Y - alpha * (grad + (epoch - 1) / (epoch + 2) * preGrad)
        preGrad = grad
        # Update matrix Q
        D = cdist(Y, Y) ** 2
        Q1 = 1 / (1 + np.log(1 + D))
        QQ1 = 1 / (1 + D)
        Q = Q1 / (np.sum(Q1) - N)
        # Compute KLD cost
        epoch = epoch + 1
        KL = -P * np.log(Q)
        KL[np.isnan(KL)] = 0
        cost = np.append(cost, np.sum(KL))
        # Update learning rate
        if (cost[-1] - cost[-2]) > 0:
            alpha = alpha * 0.99
        # Compute variation coefficient of the last three KLD costs
        if epoch > 10:
            vcc = np.var(cost[-length:]) / np.mean(cost[-length:])

    print(str(epoch - 1) + ' epochs have been computed!')
    return Y, k2