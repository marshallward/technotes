from remez import remez
import numpy as np
import matplotlib.pyplot as plt

n = 10
a, b = -0.5 * np.log(2), 0.5 * np.log(2)

func = np.exp
#p, roots, nodes, x  = remez(func, n, basis='chebyshev')
p, roots, nodes, x  = remez(func, n, domain=(a, b), basis='chebyshev')

plt.plot(x, func(x) - p(x), '+')
plt.plot(x[roots], func(x[roots]) - p(x[roots]), 'o')
plt.plot(nodes, func(nodes) - p(nodes), 'x')
plt.hlines(0, x.min(), x.max())
plt.show()
