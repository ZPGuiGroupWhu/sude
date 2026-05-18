![image](https://img.shields.io/badge/MATLAB-R2022a-brightgreen) ![image](https://img.shields.io/badge/Python-3.8-yellow) ![image](https://img.shields.io/badge/R-4.1.0-red) [![DOI](https://zenodo.org/badge/714988567.svg)](https://doi.org/10.5281/zenodo.16792257) 
# Sampling-enabled scalable manifold learning unveils the discriminative cluster structure of high-dimensional data (SUDE)
We propose a scalable manifold learning (SUDE) method that can cope with large-scale and high-dimensional data in an efficient manner. It starts by seeking a set of landmarks to construct the low-dimensional skeleton of the entire data, and then incorporates the non-landmarks into this skeleton based on the constrained locally linear embedding. This toolkit includes the main code of SUDE, and also two applications for preprocess scRNA-seq and ECG data. This paper has been published in ***Nature Machine Intelligence***, and more details can be seen https://www.nature.com/articles/s42256-025-01112-9.

![image](https://github.com/ZPGuiGroupWhu/sude/blob/master/github.png)

# How To Run
> ## MATLAB

MATLAB code of SUDE is in the ```sude_mat``` file, where the ```sude``` function provides multiple hyperparameters for user configuration as follows 
```matlab
function [Y, id_samp, para] = sude(X, varargin)
%   This function returns representation of the N by D matrix X in the lower-dimensional space and 
%   the ID of landmarks sampled by PPS. Each row in X represents an observation.
% 
%   Parameters are: 
% 
%   'NumDimensions'- A positive integer specifying the number of dimension of the representation Y. 
%                    Default: 2
%   'NumNeighbors' - A non-negative integer specifying the number of nearest neighbors for PPS to 
%                    sample landmarks. It must be smaller than N.
%                    Default: adaptive
%   'Normalize'    - Logical scalar. If true, normalize X using min-max normalization. If features in 
%                    X are on different scales, 'Normalize' should be set to true because the learning 
%                    process is based on nearest neighbors and features with large scales can override 
%                    the contribution of features with small scales. 
%                    Default: True
%   'LargeData'    - Logical scalar. If true, the data can be split into multiple blocks to avoid the problem 
%                    of memory overflow, and the gradient can be computed block by block using 'learning_l' function.                    
%                    Default: False
%   'InitMethod'   - A string specifying the method for initializing Y before manifold learning. 
%       'le'       - Laplacian eigenmaps.
%       'pca'      - Principal component analysis.
%       'mds'      - Multidimensional scaling.
%                    Default: 'le' 
%   'AggCoef'      - A positive scalar specifying the aggregation coefficient. 
%                    Default: 1.2
%   'MaxEpoch'     - Maximum number of epochs to take. 
%                    Default: 50 
```

The ```main.m``` file provides an example
```matlab
% Input data
clear;
data = csvread('benchmarks/rice.csv');
% data = textread('G:\MATLAB Drive\MATLAB\MNIST\iris.txt');

% Obtain data size and true annotations
[~, m] = size(data);
ref = data(:, m);
X = data(:, 1:m-1);
clear data

% Perform SUDE embedding
t1 = clock;
[Y, idx, para] = sude(X,'NumNeighbors',10);
t2 = clock;
disp(['Elapsed time:', num2str(etime(t2,t1)),'s']);
plotcluster2(Y, ref);
```


> ## Python

### Installation
Supported `python` versions are `3.8` and above.

This project has been uploaded to [PyPI](https://pypi.org/project/sude/), supporting direct download and installation from pypi

```
pip install sude
```

### Manual Installation

```
git clone https://github.com/ZPGuiGroupWhu/SUDE-pkg.git
cd SUDE-pkg
pip install -e .
```

The SUDE algorithm package provides the `sude` function for dimension reduction.

The description of the hyperparameters for user configuration are presented as follows

```python
def sude(
    X: np.ndarray,
    n_components: int = 2,
    *,
    n_neighbors: int = 20,
    normalize: bool = True,
    large: bool = False,
    init: Literal["le", "pca", "mds"] = "le",
    agg_coef: float = 1.2,
    max_iter: int = 50,
):
    """
    Return a lower-dimensional representation of the N by D matrix X.

    SUDE is a sampling-based scalable manifold learning method for uniform
    and discriminative embedding of large-scale and high-dimensional data. It
    first samples landmarks to construct the low-dimensional skeleton of the
    data, then incorporates non-landmark samples into this skeleton with
    constrained locally linear embedding. Each row in X represents one
    observation.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Input data matrix.
    n_components : int, default=2
        Number of dimensions in the learned embedding. Corresponds to
        ``no_dims`` in the original function interface and to the output
        dimension in the paper.
    n_neighbors : int, default=20
        Number of nearest neighbors used by PPS to sample landmarks.
        Corresponds to ``k1`` in the paper. It must be smaller than the number
        of samples when positive. Set to 0 to disable landmark sampling.
    normalize : bool, default=True
        Whether to apply min-max normalization to the input data before
        nearest-neighbor learning.
    large : bool, default=False
        Whether to use memory-bounded learning for large data.
    init : {"le", "pca", "mds"}, default="le"
        Initialization method for the embedding. Corresponds to ``initialize``
        in the original function interface and paper-style notation.
    agg_coef : float, default=1.2
        Aggregation coefficient. Corresponds to ``γ`` in the paper.
    max_iter : int, default=50
        Maximum number of optimization epochs. Corresponds to ``T_epoch`` in
        the paper.

    Returns
    -------
    Y : ndarray of shape (n_samples, n_components)
        The learned embedding.
    """
```

After installing the library, you can use the `sude` function as follows:
```python
import numpy as np
from sude import SUDE
import time
import matplotlib.pyplot as plt

# Input data
data = np.loadtxt("benchmarks/rice.csv", delimiter=",")

# Obtain data size and true annotations
m = data.shape[1]
X = data[:, :m - 1]
ref = data[:, m - 1]

# Fit a scikit-learn style estimator
start_time = time.time()
model = SUDE(
    n_components=2,
    n_neighbors=10,
    init="le",
    max_iter=50,
)
Y = model.fit_transform(X)
end_time = time.time()
print("Elapsed time:", end_time - start_time, 's')

plt.scatter(Y[:, 0], Y[:, 1], c=ref, cmap='tab10', s=4)
plt.show()
```

# Depends
> ## scRNA-seq application

argparse (≥2.0.4), assertthat (≥0.2.1), BiocGenerics (≥0.40.0), BiocSingular (≥1.10.0), ClusterR (≥1.2.5), dotCall64 (≥1.0.1), fields (≥12.5), GenomeInfoDb (≥1.30.1), GenomicRanges (≥1.46.1), geometry (≥0.4.5), ggplot2 (≥3.3.5), grid (≥4.1.0), gtools (≥3.9.2), IRanges (≥2.28.0), MatrixGenerics (≥1.6.0), mclust (≥5.4.7), parallel (≥4.1.0), prodlim (≥2019.11.13), RcppHungarian (≥0.1), readr (≥1.4.0), reshape2 (≥1.4.4), S4Vectors (≥0.30.0), scran (≥1.22.1), scuttle (≥1.4.0), Seurat (≥4.0.5), SingleCellExperiment (≥1.16.0), spam (≥2.7.0), stats4 (≥4.1.0), SummarizedExperiment (≥1.24.0), uwot (≥0.1.10)

Noted: all R packages can be installed from the [CRAN repository](https://cran.r-project.org/) or [Bioconductor](https://www.bioconductor.org/). You can also use the following R scripts to install them all.
```ruby
## Please click Tools->Global Options->Packages, change CRAN repository to a near mirror. Then, execute the following code:
## Install packages from CRAN.
install.packages(c("argparse", "assertthat", "ClusterR", "dotCall64", "fields", "geometry", "ggplot2", "gtools", "mclust", "prodlim", "RcppHungarian", "readr", "reshape2", "Seurat", "spam", "uwot"))
## Determine whether the package "BiocManager" exists, if not, install this package.
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
## Install packages from Bioconductor.
BiocManager::install(c("BiocGenerics", "BiocSingular", "GenomeInfoDb", "GenomicRanges", "IRanges", "MatrixGenerics", "S4Vectors", "scran", "scuttle", "SingleCellExperiment", "SummarizedExperiment"), force = TRUE, update = TRUE, ask = FALSE)
```

> ## ECG application

[Deep Learning Toolbox](https://ww2.mathworks.cn/products/deep-learning.html)

[Signal Processing Toolbox](https://www.mathworks.com/products/signal.html)

# Citation Request
Peng, D., Gui, Z., Wei, W. et al. Sampling-enabled scalable manifold learning unveils the discriminative cluster structure of high-dimensional data. Nat. Mach. Intell. (2025). https://doi.org/10.1038/s42256-025-01112-9
