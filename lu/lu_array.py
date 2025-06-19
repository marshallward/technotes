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
            tmp = ns.astype(U[i], U.dtype, copy=True)
            U[i] = U[pivot]
            U[pivot] = tmp

            tmp = ns.astype(P[i], P.dtype, copy=True)
            P[i] = P[pivot]
            P[pivot] = tmp

            if i > 0:
                tmp = ns.astype(L[i, :i], L.dtype, copy=True)
                L[i, :i] = L[pivot, :i]
                L[pivot, :i] = tmp

        L[i+1:, i:i+1] = U[i+1:, i:i+1] / U[i,i]

        U[i+1:, i:] = U[i+1:, i:] - L[i+1:, i:i+1] * U[i:i+1, i:]

    return P, L, U


def solve(A, b):
    """
    Solve Ax = b for 2D square matrix A and vector b.

    Parameters
    ----------
    A: array
       Square 2D floating-point array object supporting Python array API.
    b: array
       1D floating-point array object supporting the Python array API.

    Returns
    -------
    x: array
       1D floating-point array object which solves Ax = b.

    Raises
    ------
    ValueError
        If `A` is not a 2d square floating-point array or an array
        API-compatible object.
    Example
    >>> import numpy as np
    >>> from lu_array import solve
    >>> A = np.array([[2, 3], [4, 5]], dtype=float)
    >>> b = np.array([3, 3], dtype=float)
    >>> solve(A, b)
    [-3., 3.]
    """
    try:
        ns = get_namespace(A)
    except TypeError:
        raise ValueError("Matrix A is not a valid array type")

    try:
        get_namespace(b)
    except TypeError:
        raise ValueError("Vector b is not a valid array type")

    P, L, U = LU(A)

    n = A.shape[0]

    # Adjust shape of b if necessary
    # NOTE: This does not change the view of b outside of solve().
    b = ns.reshape(b, (n,))

    # Allocate the solution
    x = ns.empty(n)

    # Apply pivot
    b = P @ b

    # Forward substitution
    # NOTE: L is always diagonal!  Skip the division by L[i,i]
    # Dumb version
    for i in range(n):
        x[i] = b[i] - L[i,:i] @ x[:i]

    # Back substitution
    for i in reversed(range(n)):
        x[i] = (x[i] - U[i,i+1:] @ x[i+1:]) / U[i,i]

    return x
