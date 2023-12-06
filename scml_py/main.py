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
Y = scml(X, k1=0, large=1)
end_time = time.time()
print("Elapsed time:", end_time - start_time, 's')

plt.scatter(Y[:, 0], Y[:, 1], c=ref, cmap='tab10', s=4)
plt.show()
