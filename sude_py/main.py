import pandas as pd
import numpy as np
from sude import sude
import time
import matplotlib.pyplot as plt


# Input data
data = np.array(pd.read_csv('benchmarks/rice.csv', header=None))

# Obtain data size and true annotations
m = data.shape[1]
X = data[:, :m - 1]
ref = data[:, m - 1]

# Perform SUDE embedding
start_time = time.time()
Y = sude(X, k1=10)
end_time = time.time()
print("Elapsed time:", end_time - start_time, 's')

plt.scatter(Y[:, 0], Y[:, 1], c=ref, cmap='tab10', s=4)
plt.show()
