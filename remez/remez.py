import gmpy2
import itertools
import numpy as np
from numpy.polynomial.polynomial import Polynomial
import scipy

# What to do here?
n_grid = 4096

def remez(order, bound, func, func_hi=None):
    # The Remez algorithm iteratively generates a polynomial which minimizes
    # the error of a prescribed function.

    # Each iteration consists of two stages.
    #
    # - For n+2 prescribed points x_m:
    #
    #   - Solve for the coefficients of a Nth order polynomial p(x) whose error
    #     from a function f(x) exactly oscillates between +E and -E.
    #
    #     NOTE: This approximates - but is not necessarily - equioscillation.
    #     The errors at the xm's do oscillate, but the points are not
    #     necessarily extrema of p(x)
    #
    #   - After constructing the polynomial, do an extrema search to identify
    #     the extrema of these new points.  If not yet converged, then the
    #     number of extrema is very likely to change.
    #
    #   - Repeat the exercise until the points and the error converge.

    # Nth order polynomial in fixed doman has n+2 extrema
    n_nodes = order + 2

    # Use nodes of nth order Chebyshev polynomial as a first guess
    d_angle = 180. / n_nodes
    angles = np.linspace(0.5 * d_angle, 180. - 0.5 * d_angle, n_nodes)
    nodes = bound * scipy.special.cosdg(angles)

    # Extrema search grid
    x = np.linspace(-bound, bound, n_grid)

    for nr in range(20):
    #for nr in range(1):
    #for nr in range(2):
        #-- construct the polynomial --#

        # Construct polynomial coefficient matrix
        R = np.empty((n_nodes, n_nodes))
        R[:,:-1] = np.vander(nodes, N=n_nodes-1, increasing=True)

        # Append the error oscillation
        R[:,-1] = [(-1)**k for k in range(n_nodes)]

        # Evaluate the target function on the nodes
        f = func(nodes)

        # Solve for the polynomial coefficients and node error
        coeff = np.linalg.solve(R, f)

        # Build the polynomial
        mm_poly = Polynomial(coeff[:-1])

        # Compute error relative to high-precision exp()
        if func_hi:
            err = np.array([float(func_hi(y) - mm_poly(y)) for y in x])
        else:
            err = func(x) - mm_poly(x)

        # Alternative: Use rounded coefficients?

        ## Round each coefficient to target precision
        #rounded_coeff = np.array([float(gmpy2.next_below(gmpy2.next_above(c)))
        #        for c in coeff[:-1]])

        #mm_poly_rounded = Polynomial(rounded_coeff)

        #err = np.array([float(gmpy2.exp(y) - mm_poly_rounded(y)) for y in x])

        #-- check for convergence --#

        # Hard condition

        abs_err_bound = np.all(np.abs(err) < np.abs(coeff[-1]))

        # Soft condition: error is within 1ulp (I think?))
        # TODO: Prob not handling this correctly...
        soft_err_bound = abs(np.max(err)) - abs(coeff[-1]) < 2**-53

        # What we really want: abs(np.max(err) < 2**53)
        # but I doubt we're anywhere near that...

        # If we've converged, then we're done.
        if abs_err_bound or soft_err_bound:
            break

        #-- update nodes --#

        # Update the positions of the extrema if necessary
        # Find the roots of the error function on a dense grid.
        # TODO: Find a better method
        roots = np.where(np.sign(err[:-1]) != np.sign(err[1:]))[0]
        if roots[0] != 0:
            roots = np.insert(roots, 0, 0)
        if roots[-1] != len(err)-1:
            roots = np.append(roots, len(err)-1)

        # Now find the extrema in each segment
        # TODO: Again, write our own algorithm
        # TODO: Enforce minimum spacing if cond(R) is too high?
        nodes = np.array([
            x[l + np.abs(err[l:r+1]).argmax()]
            for l, r in itertools.pairwise(roots)
            if l != r
        ])

        n_nodes = len(nodes)

        print("len(roots) = ", len(roots))
        #print("roots = ", x[roots])
        #print("err(roots) = ", func(roots) - mm_poly(x[roots]))
        print("------")
        print("len(nodes) = ", len(nodes))
        #print("nodes = ", nodes)
        #print("err(nodes) = ", func(nodes) - mm_poly(nodes))
        print("======")
        

    #return mm_poly
    return mm_poly, roots, nodes, x
