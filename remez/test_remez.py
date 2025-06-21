from remez import remez
import numpy as np
import matplotlib.pyplot as plt

n = 12

func = np.exp
p, roots, nodes, x  = remez(n, 1, func, basis='chebyshev')
#p, roots, nodes, x  = remez(n, 1, func, basis='monomial')

plt.plot(x, func(x) - p(x), '+')
plt.plot(x[roots], func(x[roots]) - p(x[roots]), 'o')
plt.plot(nodes, func(nodes) - p(nodes), 'x')
plt.hlines(0,-1,1)
plt.show()
