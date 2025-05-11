import numpy as np
from numpy.polynomial.polynomial import Polynomial
from scipy.special import cosdg
import matplotlib.pyplot as plt
import itertools

# Construct a polynomial of order 3
order = 7

# For now, interpolate from the rescaled range [-0.5 ln(2), +0.5 ln(2)]
# TODO: Someday, look into table lookup into a smaller range.
bound = 0.5 * np.log(2.)

# Config
plot_all = False

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

# Check output
print(nodes)


#---
x = np.linspace(-bound, bound, 20000)

# Start the Remez minimax solver.
for nr in range(20):
    print(nr)

    if nr > 0:
        # Update the positions of the extrema if necessary
        # Find the roots of the error function on a dense grid.
        # TODO: Find a better method
        roots = np.where(np.sign(err[:-1]) != np.sign(err[1:]))[0]
        roots = np.insert(roots, 0, 0)
        roots = np.append(roots, len(err))
        print("roots?", roots)

        # Now find the extrema in each segment
        # TODO: Again, write our own algorithm
        # TODO: Enforce minimum spacing if cond(R) is too high?
        nodes = np.array([
            x[l + np.abs(err[l:r]).argmax()]
            for l, r in itertools.pairwise(roots)
            if l != r
        ])
        n_nodes = len(nodes)
        print("nodes?", nodes)

    # Construct the matrix
    R = np.vander(nodes, N=n_nodes-1, increasing=True)
    # TODO: Better way to do this?
    R = np.append(R, [[(-1)**k] for k in range(n_nodes)], axis=1)
    f = np.exp(nodes)
    coeff = np.linalg.solve(R, f)

    # And build the polynomial
    mm_poly = Polynomial(coeff[:-1])

    err = np.exp(x) - mm_poly(x)

    if plot_all:
        # Plot the first estimate
        plt.plot(x, err)
        plt.hlines(y=[coeff[-1], -coeff[-1]], xmin=-bound, xmax=bound)
        plt.title(f"Max Error={np.max(err)}, dE={abs(coeff[-1] - np.max(err))}")
        plt.show()

    # Hard condition
    abs_err_bound = np.all(err < coeff[-1])

    # Soft condition: error is within 1ulp (I think?)
    # TODO: Prob not handling this correctly...
    soft_err_bound = abs(np.max(err) - coeff[-1]) < 2**-53

    if abs_err_bound or soft_err_bound:
        if not plot_all:
            # Plot the first estimate
            plt.plot(x, err)
            plt.hlines(y=[coeff[-1], -coeff[-1]], xmin=-bound, xmax=bound)
            plt.title(f"Max Error={np.max(err)}, dE={abs(coeff[-1] - np.max(err))}")
            plt.show()

        break
