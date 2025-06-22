import numpy as np

def qr_decomp(A, positive_diagonals=False):
    # (This belongs in the doc file, but leave it here for now...)
    #
    # Reflect column vectors of A[i:,i:] onto basis vectors.
    #   A = H(a) R(a)

    # Brief explanation of householder reflection
    #
    # Reflection of vector x to xr across plane V (i.e. all y s.t. V.y = 0)
    #   *  xr = x + V
    #   * |x| = |x_r|
    # Then
    #   * |xr|² = |x|² + |V|² + 2 x.V and
    #   ->   |V| = -2 x.v
    # If v is the unit vector of V, then
    #   xr = x + |V| v
    #      = x - 2 (x.v) v
    #      = x - 2 v (vᵀx)
    #
    # As a transformation in matrix form:
    #   H(v) x = (I - 2 v vᵀ) x

    # Now that we have a reflection matrix H...

    # QR decomposition: since H is orthogonal:
    #   A = H0ᵀ H0 A
    #     = (H0ᵀ) (H0 A)
    # This reflects some column of A onto some basis vector e.
    # This can be used to zero out some column of A.
    #
    # Now repeat on another column
    #   (H0ᵀ H1ᵀ) (H1 H0 A)
    #   = (H0 H1) (H1 H0 A) since H is symmetric
    # First part is Q, second is R

    # TODO: reduced size support
    # Number of q basis vectors
    nq = A.shape[0]

    e = np.empty(nq)
    a = np.empty(nq)
    v = np.empty(nq)

    # Initially, A = I @ A
    Q = np.eye(nq)
    R = A.copy()

    for i in range(nq-1):
        # Construct the basis vector
        e[:] = 0.
        e[i] = 1.

        # Extract the column subvector
        a[:] = 0.
        a[i:] = R[i:,i]

        # Construct the unit plane vector
        v[:i] = 0.
        v = a + np.linalg.norm(a) * e
        v = v / np.linalg.norm(v)

        # Apply the Householder reflection to A[i,:]
        H = np.eye(nq) - 2. * np.outer(v, v)

        # Apply the latest H to Q and R
        Q = Q @ H
        R = H @ R

    if positive_diagonals:
        for i in range(min(R.shape)):
            if R[i,i] < 0:
                R[i,:] *= -1
                Q[:,i] *= -1

    return Q, R
