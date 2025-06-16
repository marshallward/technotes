import numpy as np
from numpy.polynomial.polynomial import Polynomial
from numpy.polynomial.chebyshev import chebvander, Chebyshev
from scipy.special import cosdg
import matplotlib.pyplot as plt
import itertools

import gmpy2

# Set polynomial order
order = 12

# Interpolate from the rescaled range [-0.5 ln(2), +0.5 ln(2)]
# TODO: Table lookup?
bound = 0.5 * np.log(2.)

# Config
plot_all = True
np.set_printoptions(precision=17, suppress=False)
gmpy2.get_context().precision = 113


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
x = np.linspace(-bound, bound, 10000)

e128 = [gmpy2.exp(y) for y in x]

# Start the Remez minimax solver.
for nr in range(20):

    #-- construct the polynomial --#

    # Construct the matrix
    #R = np.vander(nodes, N=n_nodes-1, increasing=True)

    scaled_nodes = nodes / bound
    R = chebvander(scaled_nodes, deg=n_nodes-2) 

    # TODO: Better way to do this?
    R = np.append(R, [[(-1)**k] for k in range(n_nodes)], axis=1)

    # Use numpy precision
    f = np.exp(nodes)

    coeff = np.linalg.solve(R, f)

    # Round each coefficient to target precision
    rounded_coeff = np.array([float(gmpy2.next_below(gmpy2.next_above(c)))
            for c in coeff[:-1]])

    # And build the polynomial
    mm_poly = Polynomial(coeff[:-1])
    mm_poly_rounded = Polynomial(rounded_coeff)
    #mm_poly_chebyshev = Chebyshev(rounded_coeff, domain=[-bound,bound])
    mm_poly_chebyshev = Chebyshev(coeff[:-1], domain=[-bound,bound])

    # Compute error relative to high-precision exp()
    #err = np.array([float(gmpy2.exp(y) - mm_poly(y)) for y in x])
    #err = np.array([float(gmpy2.exp(y) - mm_poly_rounded(y)) for y in x])
    err = np.array([float(gmpy2.exp(y) - mm_poly_chebyshev(y)) for y in x])

    #-- check for convergence --#

    # Hard condition
    abs_err_bound = np.all(err < coeff[-1])

    # Soft condition: error is within 1ulp (I think?)
    # TODO: Prob not handling this correctly...
    soft_err_bound = abs(np.max(err)) - abs(coeff[-1]) < 2**-53

    # What we really want: abs(np.max(err) < 2**53)
    # but I doubt we're anywhere near that...

    if plot_all:
        # Plot the first estimate
        plt.plot(x, err)
        plt.hlines(y=[coeff[-1], -coeff[-1]], xmin=-bound, xmax=bound)
        plt.title(f"Max Error={np.max(err)}, dE={abs(coeff[-1]) - abs(np.max(err))}")
        plt.show()
        plt.savefig(f'out_{nr}.png')

    # If we've converged, then we're done.
    if abs_err_bound or soft_err_bound:
        if not plot_all:
            # Plot the first estimate
            plt.plot(x, err)
            plt.hlines(y=[coeff[-1], -coeff[-1]], xmin=-bound, xmax=bound)
            plt.title(f"Max Error={np.max(err)}, dE={abs(coeff[-1] - np.max(err))}")
            plt.show()
            plt.savefig(f'out_{nr}.png')

        print(mm_poly)
        break

    #-- update nodes --#

    # Update the positions of the extrema if necessary
    # Find the roots of the error function on a dense grid.
    # TODO: Find a better method
    roots = np.where(np.sign(err[:-1]) != np.sign(err[1:]))[0]
    if roots[0] != 0:
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

# Show the quality of the fit?
plt.close()
#plt.plot(x, np.exp(x))
#plt.plot(x, mm_poly(x))
#plt.savefig('exp.png')
