# SUDE

SUDE is an optimized Python package for scalable manifold learning with
landmark sampling. It keeps the public API of the original `sude` package while
improving the runtime of the probability construction, gradient computation,
and non-landmark embedding steps.

The installed import package is named `sude`:

```python
from sude import SUDE, sude
```

## Installation

```bash
pip install sude
```

For best performance, install the optional Numba acceleration dependencies:
```bash
pip install "sude[accelerate]"
```
If Numba is installed, SUDE automatically uses accelerated kernels. Otherwise,
it falls back to the pure Python implementation.

## Usage

```python
import numpy as np
from sude import sude

data = np.loadtxt("benchmarks/mnist.csv", delimiter=",")
X = data[:, :-1]

Y = sude(X)
```

You can also override the original SUDE parameters:

```python
from sude import sude

Y = sude(X, no_dims=2, k1=50)
```

## API

```python
def sude(
    X,
    no_dims=2,
    k1=20,
    normalize=True,
    large=False,
    initialize="le",
    agg_coef=1.2,
    T_epoch=50,
):
    ...
```

Parameters:

- `X`: an `N` by `D` matrix. Each row represents one observation.
- `no_dims`: a positive integer specifying the dimensionality of the embedding
  `Y`. Default: `2`.
- `k1`: a non-negative integer specifying the number of nearest neighbors used
  by PPS to sample landmarks. It must be smaller than `N`. Default: `20`.
- `normalize`: whether to apply min-max normalization to `X`. If features in
  `X` are on different scales, this should usually be `True` because nearest
  neighbor learning can otherwise be dominated by large-scale features.
  Default: `True`.
- `large`: whether to use the memory-bounded learning mode, which splits
  intermediate gradient computation into blocks to reduce memory use.
  Default: `False`.
- `initialize`: the method for initializing `Y` before manifold learning.
  Supported values are `"le"` for Laplacian eigenmaps, `"pca"` for principal
  component analysis, and `"mds"` for multidimensional scaling. Default:
  `"le"`.
- `agg_coef`: a positive scalar specifying the aggregation coefficient.
  Default: `1.2`.
- `T_epoch`: maximum number of optimization epochs. Default: `50`.

The package also provides a scikit-learn compatible `SUDE` estimator. In that
API, `n_components` corresponds to `no_dims`, and `n_neighbors` corresponds to
`k1`.

## Development

Run the test suite from this directory:

```bash
python -m unittest discover -s sude/tests
```

Build the package:

```bash
python -m build
```
