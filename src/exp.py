import numpy as np
from numpy.polynomial.polynomial import Polynomial
from scipy.special import cosdg
import matplotlib.pyplot as plt

# Construct a polynomial of order 3
order = 5

# For now, interpolate from the rescaled range [-0.5 ln(2), +0.5 ln(2)]
# TODO: Someday, look into table lookup into a smaller range.
bound = 0.5 * np.log(2.)

#---

# Interpolation of an nth order polynomial requires n+1 points
n_nodes = order + 2

# Initialize with Chebyshev nodes, which should approximate the extrema and
# capture oscillation of the error.
# First and last point are extrema, separated by 0.5 dTheta, but the n zeros
# are separated by dTheta.
# 180 = 0.5 dTheta + n dTheta + 0.5 dTheta

d_angle = 180. / n_nodes

angles = np.linspace(0.5 * d_angle, 180. - 0.5 * d_angle, n_nodes)
nodes = bound * cosdg(angles)

## Check output
#print(angles)
#print(nodes)


#---

# Start the Remez minimax solver.  Set up the first iteration.

R = np.vander(nodes, N=n_nodes-1, increasing=True)
# TODO: Better way to do this?
R = np.append(R,[[(-1)**k] for k in range(n_nodes)],axis=1)

f = np.exp(nodes)

# Now solve for the coefficients
coeff = np.linalg.solve(R, f)

## Check output
#print(R)
#print(f)
#print(coeff)

#---

# Check results

# Plot the first estimate
x = np.linspace(-bound, bound, 100)
mm_poly = Polynomial(coeff[:-1])

plt.plot(x, np.exp(x) - mm_poly(x))
plt.show()

#---

# Try a second iteration of Remez

# 1. Find the new extrema






