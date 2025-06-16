import numpy as np
from numpy.polynomial.polynomial import Polynomial
from numpy.polynomial.chebyshev import chebvander, Chebyshev
from scipy.special import cosdg
import matplotlib.pyplot as plt
import itertools

import gmpy2
import mpmath

# Set polynomial order
order = 12

# Interpolate from the rescaled range [-0.5 ln(2), +0.5 ln(2)]
# TODO: Table lookup?
bound = 0.5 * np.log(2.)
#bound = 1

# Config
plot_all = True
np.set_printoptions(precision=17, suppress=False)
gmpy2.get_context().precision = 113
mpmath.mp.dps = 34  # r128?

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
print("initial nodes:")
print(nodes)
print("start remez")

#---

def chebyshev_eval(coeffs, x_scaled):
    # coeffs are [c0, c1, c2, ...]
    T0 = gmpy2.mpfr(1)
    T1 = x_scaled
    result = coeffs[0] * T0
    if len(coeffs) > 1:
        result += coeffs[1] * T1
    for k in range(2, len(coeffs)):
        Tk = 2 * x_scaled * T1 - T0
        result += coeffs[k] * Tk
        T0, T1 = T1, Tk
    return result


#---
x = np.linspace(-bound, bound, 10000)

# Start the Remez minimax solver.
for nr in range(20):

    #-- construct the polynomial --#

    # Construct the matrix
    #R = np.vander(nodes, N=n_nodes-1, increasing=True)

    #scaled_nodes = nodes / bound

    #R = chebvander(scaled_nodes, deg=n_nodes-2) 
    ## TODO: Better way to do this?
    #R = np.append(R, [[(-1)**k] for k in range(n_nodes)], axis=1)

    scaled_nodes = [gmpy2.mpfr(n / bound) for n in nodes]

    R_mp = []
    for xi in scaled_nodes:
        row = [gmpy2.cos(k * gmpy2.acos(xi)) for k in range(n_nodes-1)]
        row.append(gmpy2.mpfr((-1)**len(R_mp)))
        R_mp.append(row)

    # Use numpy precision
    #f = np.exp(nodes)
    f_mp = [gmpy2.exp(n) for n in nodes]

    #coeff = np.linalg.solve(R, f)
    ## Round each coefficient to target precision
    #rounded_coeff = np.array([float(gmpy2.next_below(gmpy2.next_above(c)))
    #        for c in coeff[:-1]])

    R_mat = mpmath.matrix(R_mp)
    f_vec = mpmath.matrix(f_mp)
    coeff_mp = mpmath.lu_solve(R_mat, f_vec)

    ## And build the polynomial
    #mm_poly = Polynomial(coeff[:-1])
    #mm_poly_rounded = Polynomial(rounded_coeff)
    ##mm_poly_chebyshev = Chebyshev(rounded_coeff, domain=[-bound,bound])
    #mm_poly_chebyshev = Chebyshev(coeff[:-1], domain=[-bound,bound])

    # Compute error relative to high-precision exp()
    #err = np.array([float(gmpy2.exp(y) - mm_poly(y)) for y in x])
    #err = np.array([float(gmpy2.exp(y) - mm_poly_rounded(y)) for y in x])
    #err = np.array([float(gmpy2.exp(y) - mm_poly_chebyshev(y)) for y in x])

    x_mp = [gmpy2.mpfr(xx) for xx in x]  # scaled x
    err_mp = [float(gmpy2.exp(xx) - chebyshev_eval(coeff_mp[:-1], xx/bound)) for xx in x_mp]

    #-- check for convergence --#
    err = np.array([abs(float(e)) for e in err_mp])
    print("pointwise error:")
    print(err)

    print("remez error:")
    print(coeff_mp)
    print(coeff_mp[1], coeff_mp[len(coeff_mp)-1], coeff_mp[-1])

    # Hard condition
    abs_err_bound = np.all(np.array([e < float(coeff_mp[len(coeff_mp)-1]) for e in err]))

    # Soft condition: error is within 1ulp (I think?)
    # TODO: Prob not handling this correctly...
    #soft_err_bound = abs(np.max(err)) - abs(float(coeff_mp[-1])) < 2**-52
    soft_err_bound = np.all(np.array([e - float(coeff_mp[len(coeff_mp)-1]) < 2**-52 for e in err]))

    print(abs_err_bound or soft_err_bound)

    # What we really want: abs(np.max(err) < 2**53)
    # but I doubt we're anywhere near that...

    #if plot_all:
    #    # Plot the first estimate
    #    plt.plot(x, err)
    #    plt.hlines(y=[float(coeff_mp[-1]), -float(coeff_mp[-1])], xmin=-bound, xmax=bound)
    #    plt.title(f"Max Error={np.max(err)}, dE={abs(float(coeff_mp[-1])) - abs(np.max(err))}")
    #    plt.show()
    #    plt.savefig(f'out_{nr}.png')

    # If we've converged, then we're done.
    if abs_err_bound or soft_err_bound:
        #if not plot_all:
        #    # Plot the first estimate
        #    plt.plot(x, err)
        #    plt.hlines(y=[float(coeff_mp[-1]), -float(coeff_mp[-1])], xmin=-bound, xmax=bound)
        #    plt.title(f"Max Error={np.max(err)}, dE={abs(float(coeff_mp[-1]) - np.max(err))}")
        #    plt.show()
        #    plt.savefig(f'out_{nr}.png')

        #print(mm_poly)
        break

    #-- update nodes --#

    ## Update the positions of the extrema if necessary
    ## Find the roots of the error function on a dense grid.
    ## TODO: Find a better method
    #err = np.array([float(e) for e in err_mp])
    #roots = np.where(np.sign(err[:-1]) != np.sign(err[1:]))[0]
    #if roots[0] != 0:
    #    roots = np.insert(roots, 0, 0)
    #roots = np.append(roots, len(err))
    #print("roots?", roots)

    ## Now find the extrema in each segment
    ## TODO: Again, write our own algorithm
    ## TODO: Enforce minimum spacing if cond(R) is too high?
    #nodes = np.array([
    #    x[l + np.abs(err[l:r]).argmax()]
    #    for l, r in itertools.pairwise(roots)
    #    if l != r
    #])

    # ChatGPT is our friend
    # Compute sign of mpfr errors robustly
    err_sign = np.array([gmpy2.cmp(e, gmpy2.mpfr(0)) for e in err_mp])
    
    # Find sign changes
    roots = np.where(err_sign[:-1] != err_sign[1:])[0]
    
    # Ensure endpoints included
    if roots[0] != 0:
        roots = np.insert(roots, 0, 0)
    roots = np.append(roots, len(err_sign))
    
    print("roots?", roots)
    
    # Now find extrema in each segment (still safe to use float here)
    # because we're only using argmax() of abs(error)
    err_abs = np.array([float(abs(e)) for e in err_mp])
    
    nodes = np.array([
        x[l + np.argmax(err_abs[l:r])]
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
