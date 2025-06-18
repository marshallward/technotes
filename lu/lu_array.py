from array_api_compat import get_namespace

def LU(A):
    try:
        ns = get_namespace(A)
    except TypeError:
        raise ValueError("Input is not a valid array type")

    A_shape = ns.shape(A)
    if len(A_shape) != 2 or A_shape[0] != A_shape[1]:
        raise ValueError("Input must be a square 2D array")

    n = A_shape[0]

    U = ns.copy(A)
    L = ns.eye(n)
    P = ns.eye(n)

    for i in range(n-1):
        # Partial pivot to row by maximum diagonal

        col = ns.abs(U[i:,i])
        offset = ns.argmax(col)
        pivot = i + offset

        # Swap pivot row
        if pivot != i and not ns.isclose(col[offset], col[0], atol=1e-12):

            # NOTE: Array API does not require U[[i,pivot]] = U[[pivot, i]]
            tmp = ns.copy(U[i])
            U[i] = U[pivot]
            U[pivot] = tmp
            
            tmp = ns.copy(P[i])
            P[i] = P[pivot]
            P[pivot] = tmp

            if i > 0:
                tmp = ns.copy(L[i, :i])
                L[i, :i] = L[pivot, :i]
                L[pivot, :i] = tmp

        L[i+1:, i:i+1] = U[i+1:, i:i+1] / U[i,i]

        U[i+1:, i:] = U[i+1:, i:] - L[i+1:, i:i+1] * U[i:i+1, i:]

    return P, L, U
