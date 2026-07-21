import gmpy2
import itertools
import numpy as np
from numpy.polynomial.polynomial import Polynomial
from numpy.polynomial.chebyshev import chebvander, Chebyshev
import scipy

# Precision used for high-precision target and polynomial evaluation.
gmpy2.get_context().precision = 200

# What to do here?
n_grid = 16385
n_iter = 20


# NOTE: This is ChatGPT brainfart, can probably be integrated into roots[:]
# NOTE: If too high (say =10) then it hurts convergence.
def suppress_nearby_roots(roots, min_spacing=2):
    if len(roots) == 0:
        return roots
    filtered = [roots[0]]
    for r in roots[1:]:
        if r - filtered[-1] >= min_spacing:
            filtered.append(r)
    return np.array(filtered)


def polyval_hi(coeff, x):
    """Evaluate a monomial-basis polynomial in gmpy2 precision."""
    x = gmpy2.mpfr(x)
    value = gmpy2.mpfr(0)

    for c in reversed(coeff):
        value = value * x + gmpy2.mpfr(float(c))

    return value


def chebval_hi(coeff, x_scaled):
    """Evaluate a Chebyshev series in gmpy2 precision using Clenshaw."""
    x_scaled = gmpy2.mpfr(x_scaled)
    b_kplus1 = gmpy2.mpfr(0)
    b_kplus2 = gmpy2.mpfr(0)

    for c in reversed(coeff[1:]):
        b_k = 2 * x_scaled * b_kplus1 - b_kplus2 + gmpy2.mpfr(float(c))
        b_kplus2 = b_kplus1
        b_kplus1 = b_k

    return x_scaled * b_kplus1 - b_kplus2 + gmpy2.mpfr(float(coeff[0]))

def remez(func, order, domain=None, func_hi=None, basis='monomial', debug=False):
    # The Remez algorithm iteratively generates a polynomial which minimizes
    # the error of a prescribed function.

    if basis not in ('monomial', 'chebyshev'):
        raise ValueError('Unknown basis: {basis!r}')

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
    expected_nodes = n_nodes

    # Use nodes of nth order Chebyshev polynomial as a first guess
    d_angle = 180. / n_nodes
    angles = np.linspace(0.5 * d_angle, 180. - 0.5 * d_angle, n_nodes)

    # Rescale Chebyshev nodes from [-1,1] to [a,b]
    scaled_nodes = scipy.special.cosdg(angles)
    if domain:
        median = 0.5 * (domain[0] + domain[1])
        nodes = median + 0.5 * (domain[1] - domain[0]) * scaled_nodes
    else:
        nodes = scaled_nodes

    # Extrema search grid
    if domain:
        x = np.linspace(domain[0], domain[1], n_grid)
    else:
        x = np.linspace(-1, 1, n_grid)

    for ni in range(n_iter):
        if debug:
            print(ni)

        #-- construct the polynomial --#

        # Construct polynomial coefficient matrix
        R = np.empty((n_nodes, n_nodes))

        if basis == 'monomial':
            R[:,:-1] = np.vander(nodes, N=n_nodes-1, increasing=True)
        elif basis == 'chebyshev':
            R[:,:-1] = chebvander(scaled_nodes, deg=n_nodes-2)

        # Append the error oscillation
        R[:,-1] = [(-1)**k for k in range(n_nodes)]

        # Evaluate the target function on the nodes
        f = func(nodes)

        # Solve for the polynomial coefficients and node error
        coeff = np.linalg.solve(R, f)

        # Build the polynomial
        if basis == 'monomial':
            mm_poly = Polynomial(coeff[:-1])
        elif basis == 'chebyshev':
            mm_poly = Chebyshev(coeff[:-1], domain=domain)

        # Compute error relative to high-precision exp()
        if func_hi:
            poly_coeff = coeff[:-1]

            if basis == 'monomial':
                err = np.array([
                    float(func_hi(y) - polyval_hi(poly_coeff, y))
                    for y in x
                ])
            elif basis == 'chebyshev':
                x_scaled = (
                    2. * x / (domain[1] - domain[0])
                    - (domain[1] + domain[0]) / (domain[1] - domain[0])
                    if domain is not None else x
                )
                err = np.array([
                    float(func_hi(y) - chebval_hi(poly_coeff, ys))
                    for y, ys in zip(x, x_scaled)
                ])
        else:
            err = func(x) - mm_poly(x)

        # Alternative: Use rounded coefficients?

        ## Round each coefficient to target precision
        #rounded_coeff = np.array([float(gmpy2.next_below(gmpy2.next_above(c)))
        #        for c in coeff[:-1]])

        #mm_poly_rounded = Polynomial(rounded_coeff)

        #err = np.array([float(gmpy2.exp(y) - mm_poly_rounded(y)) for y in x])

        #-- update nodes --#

        # Update the positions of the extrema if necessary

        # Find the roots of the error function on a dense grid.

        #--#

        # First make exact zeros slightly positive
        # (Actually, GPT says to try preserve adjacent sign...)

        err_adj = err.copy()
        zero = err_adj == 0.0

        # Propagate the previous nonzero value forward.
        for i in range(1, len(err_adj)):
            if zero[i]:
                err_adj[i] = err_adj[i - 1]

        # Handle any leading zeros by propagating backward.
        for i in range(len(err_adj) - 2, -1, -1):
            if err_adj[i] == 0.0:
                err_adj[i] = err_adj[i + 1]

        roots = np.where(np.sign(err_adj[:-1]) != np.sign(err_adj[1:]))[0]

        # Does this work?
        roots = suppress_nearby_roots(roots)

        if roots.size == 0:
            # No roots!  This could be good or bad...
            # For now, assume its FPE noise (i.e. "good news")
            break

        # Append the endpoints, to define endpoint intervals
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

        if domain:
            a, b = domain
            scaled_nodes = 2. * nodes / (b - a) - (b + a) / (b - a)
        else:
            scaled_nodes = nodes

        #-- check for convergence --#

        # Hard condition

        abs_err_bound = np.all(np.abs(err) < np.abs(coeff[-1]))

        # Soft condition: error is within 1ulp (I think?))
        # TODO: Prob not handling this correctly...
        soft_err_bound = abs(np.max(err)) - abs(coeff[-1]) <= 2**-51

        # What we really want: abs(np.max(err) < 2**53)
        # but I doubt we're anywhere near that...

        # If we've converged, then we're done.
        if abs_err_bound or soft_err_bound:
            break

        # I used to just change the order:
        #n_nodes = len(nodes)

        # ... but GPT said this is dangerous.  Ok?
        if len(nodes) != expected_nodes:
            print(
                f"iter {ni}: "
                f"roots={len(roots)} "
                f"nodes={len(nodes)} "
                f"expected={order+2}"
            )

            # But break if you just want the coefs
            raise RuntimeError(
                f"Remez exchange failed at iteration {ni}: "
                f"expected {expected_nodes} extrema, found {len(nodes)}"
            )

        if debug:
            print("len(roots) = ", len(roots))
            #print("roots = ", x[roots])
            #print("err(roots) = ", func(roots) - mm_poly(x[roots]))
            print("------")
            print("len(nodes) = ", len(nodes))
            #print("nodes = ", nodes)
            #print("err(nodes) = ", func(nodes) - mm_poly(nodes))
            print("======")


    if debug:
        return mm_poly, roots, nodes, x
    else:
        return mm_poly
