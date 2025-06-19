import numpy as np
import scipy

from lu_array import solve

A = np.array([[2., 3.], [4., 5.]])
b = np.array([3., 3.])

x = solve(A,b)
print(x)


A = scipy.linalg.hilbert(8)
b = np.linspace(0., 1., 8)

x = solve(A,b)
print(x)

# ref
print(scipy.linalg.solve(A,b))
