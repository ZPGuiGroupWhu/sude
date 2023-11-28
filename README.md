![image](https://img.shields.io/badge/MATLAB-R2022a-brightgreen) ![image](https://img.shields.io/badge/Python-3.11-blue) ![image](https://img.shields.io/badge/R-4.1.0-red) 
# scML: scalable manifold learning by uniform landmark sampling and constrained locally linear embedding
We propose a scalable manifold learning (scML) method that can cope with large-scale and high-dimensional data in an efficient manner. It starts by seeking a set of landmarks to construct the low-dimensional skeleton of the entire data, and then incorporates the non-landmarks into this skeleton based on the constrained locally linear embedding. This toolkit includes the main code of scML, and also two applications for preprocess scRNA-seq and ECG data.

![image](https://github.com/ZPGuiGroupWhu/scml/blob/main/github.png)

# How To Run
> **MATLAB**

The MATLAB code of scML is in the 'scml' file, where the 'scml' function provides numerous hyperparameters for user configuration as follows 
```matlab
function [Y, id_samp, para] = scml(X, varargin)
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
%   'TolVcc'       - Termination tolerance for variation coefficient of the last three KLD costs. 
%                    Default: 1e-7    
```

Open the 'main.m' file, we provide an example
```matlab
% Input data
clear;
data = csvread('benchmarks/ds3.csv');

% Obtain data size and true annotations
[~, m] = size(data);
ref = data(:, m);
X = data(:, 1:m-1);
clear data

% Perform scML embedding
t1 = clock;
[Y, idx, para] = scml(X);
t2 = clock;
disp(['Elapsed time:', num2str(etime(t2,t1)),'s']);
plotcluster2(Y, ref);
```


> **Python**

The Python code of scML is in the 'scml_py' file, where the 'scml' function provides numerous hyperparameters for user configuration as follows
```python
def scml(
    X,
    no_dims = 2,
    k1 = 20,
    normalize = True,
    large = False,
    initialize = 'le',
    agg_coef = 1.2,
    T_epoch = 50,
    T_vcc = 1e-7
):
"""
    This function returns representation of the N by D matrix X in the lower-dimensional space. Each row in X
    represents an observation.

    Parameters are:

    'no_dims'      - A positive integer specifying the number of dimension of the representation Y.
                   Default: 2
    'k1'           - A non-negative integer specifying the number of nearest neighbors for PPS to
                   sample landmarks. It must be smaller than N.
                   Default: adaptive
    'normalize'    - Logical scalar. If true, normalize X using min-max normalization. If features in
                   X are on different scales, 'Normalize' should be set to true because the learning
                   process is based on nearest neighbors and features with large scales can override
                   the contribution of features with small scales.
                   Default: True
    'large'        - Logical scalar. If true, the data can be split into multiple blocks to avoid the problem
                   of memory overflow, and the gradient can be computed block by block using 'learning_l' function.
                   Default: False
    'initialize'   - A string specifying the method for initializing Y before manifold learning.
        'le'       - Laplacian eigenmaps.
        'pca'      - Principal component analysis.
        'mds'      - Multidimensional scaling.
                   Default: 'le'
    'agg_coef'     - A positive scalar specifying the aggregation coefficient.
                   Default: 1.2
    'T_epoch'      - Maximum number of epochs to take.
                   Default: 50
    'T_vcc'        - Termination tolerance for variation coefficient of the last three KLD costs.
                   Default: 1e-7
"""
```

Open the 'main.py' file, we provide an example
```python
import pandas as pd
import numpy as np
from scml import scml
import time
import matplotlib.pyplot as plt

# Input data
data = np.array(pd.read_csv('benchmarks/ds3.csv', header=None))

# Obtain data size and true annotations
m = data.shape[1]
X = data[:, :m - 1]
ref = data[:, m - 1]

# Perform scML embedding
start_time = time.time()
Y = scml(X)
end_time = time.time()
print("Elapsed time:", end_time - start_time, 's')

plt.scatter(Y[:, 0], Y[:, 1], c=ref, cmap='tab10', s=4)
plt.show()
```

# Depends
> **scRNA-seq application**

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

> **ECG application**

[Deep Learning Toolbox](https://ww2.mathworks.cn/products/deep-learning.html)

[Signal Processing Toolbox](https://www.mathworks.com/products/signal.html)
