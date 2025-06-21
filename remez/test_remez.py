from remez import remez
import numpy as np
import matplotlib.pyplot as plt

p, roots, nodes, x  = remez(10, 1, np.sin)

plt.plot(x, np.sin(x) - p(x))
plt.plot(x[roots], np.sin(x[roots]) - p(x[roots]), 'o')
plt.plot(nodes, np.sin(nodes) - p(nodes), 'x')
plt.hlines(0,-1,1)
plt.show()
