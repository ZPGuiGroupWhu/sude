The Python versions of SUDE can be found at: https://github.com/ZPGuiGroupWhu/SUDE-pkg.

## Project layout

The project now follows the structure of the
`scikit-learn-contrib/project-template`:

```text
.
|-- .github/workflows/
|-- benchmarks/
|-- doc/
|-- examples/
|-- image/
|-- sude/
|   |-- __init__.py
|   |-- _learning_utils.py
|   |-- _numba_kernels.py
|   |-- _sude.py
|   |-- _version.py
|   `-- learning.py
|-- tests/
|-- pyproject.toml
`-- README.md
```

## Installation
Supported `python` versions are `3.8` and above.

This project has been uploaded to [PyPI](https://pypi.org/project/sude/), supporting direct download and installation from pypi

```
pip install sude
```

Numba-accelerated kernels are installed by default. SUDE enables them
automatically when both the unique input sample count and the landmark count
are large enough.

The default thresholds are:

```python
NUMBA_AUTO_MIN_SAMPLES = 3000
NUMBA_AUTO_MIN_LANDMARKS = 512
```

You can adjust them before fitting:

```python
import sude.learning as sude_learning

sude_learning.NUMBA_AUTO_MIN_SAMPLES = 5000
sude_learning.NUMBA_AUTO_MIN_LANDMARKS = 1024
```

Both values must be positive integers. If either value is invalid, SUDE falls
back to using numba whenever numba is installed.

### Manual installation

```
git clone https://github.com/ZPGuiGroupWhu/SUDE-pkg.git
cd SUDE-pkg
pip install -e .
```

## How to run

The package now exposes both a scikit-learn style estimator class and a
function wrapper with matching parameter names.

### Estimator interface

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

The estimator provides the familiar API:

```python
model = SUDE(n_components=2, n_neighbors=10, init="le")
Y = model.fit_transform(X)
Y_new = model.transform(X_new)
```

### Function interface

The function entry point uses the same sklearn-style parameter names as the
estimator:

```python
from sude import sude

Y = sude(X, n_components=2, n_neighbors=10, init="le", max_iter=50)
```

For readers comparing with the paper or original function interface,
``n_components`` corresponds to ``no_dims``, ``n_neighbors`` corresponds to
``k1``, ``init`` corresponds to ``initialize``, and ``max_iter`` corresponds to
``T_epoch``.

Run the packaged example with:

```bash
uv run python examples/plot_sude_embedding.py
```

Run the test suite with:

```bash
uv run python -m unittest discover -s tests
```

## Citation request
Peng, D., Gui, Z., Wei, W. et al. Sampling-enabled scalable manifold learning unveils the discriminative cluster structure of high-dimensional data. Nat. Mach. Intell. (2025). https://doi.org/10.1038/s42256-025-01112-9
