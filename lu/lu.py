import numpy as np
import scipy

### Somewhat arbitrary example
#np.random.seed(0)
#A = np.random.rand(4,4)

## Textbook example
#A = np.array([[2., 1., -1.], [-3., -1., 2.], [-2., 1., 2.]])

# A high condition number example
#A = scipy.linalg.hilbert(8)
A = scipy.linalg.hilbert(4)

## A pivoting example
#A = np.array([[2., 1., -1.], [4., 2., -2.], [-2., 1., 2.]])

print("A = ")
print(A)

assert(A.shape[0] == A.shape[1])
n = A.shape[0]

E = np.empty((n,n))

# To compute the echelon form of A, we march across the columns of A.
# For diagonal point A[i,i]:
#   For row j >= i:
#       E[i:,j] = A[i:,j] - (a[i,:] / a[i,i]) * A[i:,i]
#
#   This sets column A[i:,i+1:] to zero, and modifies the rest.

# Compute the echelon (upper triangular) form of A

# This is like the "first term" in each element of E
E = A.copy()
P = np.eye(n)
L = np.eye(n)

for i in range(n-1):
    #if debug:
    #    print("i=", i)
    #    # Collect the i column coefficient, weighted by the i diagonal
    #    c = E[i+1:,i:i+1] / E[i,i]
    #    print("c=")
    #    print(c)

    #    # Compute the diff row
    #    r = E[i:i+1,i:]
    #    print("r=")
    #    print(r)
    #
    #    # Apply to all i+1 rows, which should force E[i+1:,i] to zero
    #    E[i+1:,i:] = E[i+1:,i:] - c * r
    #    print("E=")
    #    print(E)

    ## Probably faster to do as a single line
    #E[i+1:,i:] = E[i+1:,i:] - (E[i+1:,i:i+1] / E[i,i]) * E[i:i+1,i:]

    # Partial pivot
    # - Search final i rows (including i) for the max coefficient.
    #   - Not that we include the offset to global index
    # - Swap the rows, so that we divide by this largest value (?)
    #pivot = np.argmax(E[i:,i]) + i
    #if pivot != i:

    col = np.abs(E[i:,i])
    pivot_offset = np.argmax(col)
    pivot = i + pivot_offset

    if pivot != i and not np.isclose(col[pivot_offset], col[0], atol=1e-12):
        E[[i, pivot]] = E[[pivot, i]]
        P[[i, pivot]] = P[[pivot, i]]

        # swap the previous L rows up to column i
        # huh??
        if i > 0:
            L[[i, pivot], :i] = L[[pivot, i], :i]
        #L[[i, pivot], :i] = L[[pivot, i], :i]

    ## This expression also loops over the E[:i,:].
    ## E[:i,:] = 0, so it should be numerically equivalent.  It may be faster,
    ## despite the extra work, due to pipelining and vectorization.
    ## But I haven't yet tested.
    #E[i+1:,:] = E[i+1:,:] - (E[i+1:,i:i+1] / E[i,i]) * E[i:i+1,:]

    ## This equivalent but relies on implict reshape of E[i,:] from (n,) to (1,n)
    ## (I don't think it works if I do pivot swaps)
    #E[i+1:,:] = E[i+1:,:] - (E[i+1:,i:i+1] / E[i,i]) * E[i,:]

    # This method moves the correction terms into a separate matrix (L)
    L[i+1:,i:i+1] = E[i+1:,i:i+1] / E[i,i]

    # Then independently apply it to E
    E[i+1:,:] = E[i+1:,:] - L[i+1:,i:i+1] * E[i:i+1,:]


# E is in echelon form
print("U = ")
print(E)
print("L = ")
print(L)
print("P = ")
print(P)


# Check the answer

if np.allclose(P @ L @ E, A):
    print("P L U = A")
else:
    print("ERROR! P L U /= A!")
    print("A=")
    print(A)
    print("P @ L @ U =")
    print(P @ L @ E)
