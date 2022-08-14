
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import pandas as pd


df = pd.read_csv('dist_puck.dat', header=None, delimiter="\s+")
X = df.columns.values
# print(X)
Y = df.index.values
# print(Y)
Z = df.values
# print(Z)
x, y = np.meshgrid(X, Y)
plt.contourf(x, y, Z)
# plt.show()

fig, ax = plt.subplots()
CS = ax.contourf(x, y, Z, cmap="RdBu")
ax.set_aspect('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')

fig.colorbar(CS, format="%.3f")
plt.show()
