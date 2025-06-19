from array_api_compat import get_namespace

def LU(A):
    """
    Perform LU decomposition with partial pivoting.

    Parameters
    ----------
    A: array
       Square 2D floating-point array object supporting the Python Array API.

    Returns
    -------
    P: array
       Permutation matrix such that P @ A = L @ U.
    L: array
       Lower triangular matrix with unit diagonal.
    U: array
       Upper triangular matrix.

    Raises
    ------
    ValueError
        If `A` is not a 2d square floating-point array or an array
        API-compatible object.

    Example
    >>> import numpy as np
    >>> from lu_array import LU
    >>> A = np.array([[2, 3], [4, 5]], dtype=float)
    >>> P, L, U = LU(A)
    >>> np.allclose(P @ A, L @ U)
    True
    """
    try:
        ns = get_namespace(A)
    except TypeError:
        raise ValueError("Input is not a valid array type")

    if not str(A.dtype).startswith("float"):
        A = ns.astype(A, ns.float64)

    if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("Input must be a square 2D array")

    n = A.shape[0]

    # ns.copy(A) is valid, but not yet supported in PyTorch
    U = ns.astype(A, A.dtype, copy=True)
    L = ns.eye(n, dtype=A.dtype)
    P = ns.eye(n, dtype=A.dtype)

    for i in range(n-1):
        # Partial pivot to row by maximum diagonal

        col = ns.abs(U[i:,i])
        offset = ns.argmax(col)
        pivot = i + offset

        # Swap pivot row
        if pivot != i and abs(float(col[offset])) - abs(float(col[0])) > 1e-12:

            # NOTE: Array API does not require U[[i,pivot]] = U[[pivot, i]]
            #   - JAX cannot even do U[i] = U[pivot]!  Immutable!!
            tmp = ns.astype(U[i], U.dtype, copy=True)
            U = U.at[i].set(U[pivot])
            U = U.at[pivot].set(tmp)

            tmp = ns.astype(P[i], P.dtype, copy=True)
            P = P.at[i].set(P[pivot])
            P = P.at[pivot].set(tmp)

            # i = 0 does nothing but it may be a slight optimization to skip
            if i > 0:
                tmp = ns.astype(L[i, :i], L.dtype, copy=True)
                L = L.at[i, :i].set(L[pivot, :i])
                L = L.at[pivot, :i].set(tmp)

        #L[i+1:, i:i+1] = U[i+1:, i:i+1] / U[i,i]
        L = L.at[i+1:, i:i+1].set(U[i+1:, i:i+1] / U[i,i])

        #U[i+1:, i:] = U[i+1:, i:] - L[i+1:, i:i+1] * U[i:i+1, i:]
        U = U.at[i+1:, i:].set(U[i+1:, i:] - L[i+1:, i:i+1] * U[i:i+1, i:])

    return P, L, U
