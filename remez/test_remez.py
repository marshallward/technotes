from remez import remez
import numpy as np
import matplotlib.pyplot as plt
import gmpy2

n = 9
a, b = -0.5 * np.log(2), 0.5 * np.log(2)

def g(x):
    # Use masked division
    return np.divide(np.expm1(x), x,
        out=np.ones_like(x),
        where=(x != 0.0),
    )

    return result

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
#p = remez(g, n, domain=(a,b), basis='chebyshev')
p = remez(g, n, domain=(a,b), basis='monomial')

x = np.linspace(a, b, 1024)
#plt.plot(x, p(x) - func(x))
plt.plot(x, x * p(x) - x * g(x))
plt.show()

# Print coefficients
print(f"real(real64), parameter :: c(0:{len(p.coef)-1}) = [ &")

for i, c in enumerate(p.coef):
    end = ", &" if i < len(p.coef)-1 else " ]"
    print(f"    {c:.17e}_real64{end}")
