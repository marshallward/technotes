import numpy as np

def qr_decomp(A, mode='full', positive_diagonals=False):
    """
    Compute the QR decomposition of a matrix A using Householder reflections.

    Parameters
    ----------
    A : (m, n) array_like
        A real-valued 2D array to decompose.
    mode: {'full', 'reduced'}, optional
        Return full matrices, even if problem is overdetermined.
        'full': Q is mxm, R is mxn with zero-rows retained.
        'reduced': Q is mxn, R is nxn with zero-rows trimmed.
    positive_diagonals : bool, optional
        If True, ensures that all diagonal elements of R are positive.

    Returns
    -------
    Q : (m, m) ndarray
        An orthogonal matrix.
    R : (m, n) ndarray
        An upper triangular matrix such that Q @ R = A.
    """
    # TODO: reduced size support
    # Number of q basis vectors
    m, n = A.shape

    # m = size of column vectors
    # n = number of column vectors

    e = np.empty(m)
    a = np.empty(m)
    v = np.empty(m)

    # Initially, A = I @ A
    Q = np.eye(m)
    R = A.copy()

    for i in range(min(m,n)):
        # Construct the basis vector
        # NOTE: For some reason it's faster to use np.zeros()
        e = np.zeros(m)
        e[i] = 1.

        # Extract the column subvector
        a[:i] = 0.
        a[i:] = R[i:,i]

        a_norm = np.linalg.norm(a)
        if a_norm == 0:
            # Skip zero-norm columns
            # I guess they're already zero?
            continue

        # Construct the unit plane vector
        v = a + np.sign(a[i]) * a_norm * e
        v[:] = v / np.linalg.norm(v)

        # Apply the Householder reflection to A[i,:]
        H = np.eye(m) - 2. * np.outer(v, v)

        # Apply the latest H to Q and R
        Q = Q @ H
        R = H @ R

    if positive_diagonals:
        for i in range(min(R.shape)):
            if R[i,i] < 0:
                R[i,:] *= -1
                Q[:,i] *= -1

    if mode == 'reduced' and m > n:
        return Q[:,:n], R[:n, :]
    else:
        return Q, R
