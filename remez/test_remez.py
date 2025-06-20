from remez import remez
import numpy as np
import matplotlib.pyplot as plt

p, roots, nodes, x  = remez(12, 1, np.sin)

plt.plot(x, np.sin(x) - p(x))
plt.plot(nodes, np.sin(nodes) - p(nodes), 'x')
plt.show()
