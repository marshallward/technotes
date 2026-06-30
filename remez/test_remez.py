from remez import remez
import numpy as np
import matplotlib.pyplot as plt
import gmpy2

n = 9
a, b = -0.5 * np.log(2), 0.5 * np.log(2)

def g(x):
    return np.where(x == 0, 1.0, np.expm1(x) / x)

def g_hi(x):
    return gmpy2.expm1(x) / x if x != 0 else gmpy2.mpfr(1)

#func = np.sin
#func = np.exp
#p, roots, nodes, x  = remez(func, n, domain=(a, b), basis='chebyshev')
##p, roots, nodes, x  = remez(func, n, domain=(a, b), basis='monomial')
#
#plt.plot(x, func(x) - p(x), '+')
#plt.plot(x[roots], func(x[roots]) - p(x[roots]), 'o')
#plt.plot(nodes, func(nodes) - p(nodes), 'x')
#plt.hlines(0, x.min(), x.max())
#plt.show()

#p = remez(func, n, domain=(a,b), basis='chebyshev')
#p = remez(g, 10, domain=(a,b), basis='chebyshev')
p = remez(g, 11, domain=(a,b), basis='chebyshev')

x = np.linspace(a, b, 1024)
#plt.plot(x, p(x) - func(x))
plt.plot(x, x * p(x) - x * g(x))
plt.show()
