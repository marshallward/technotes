import numpy as np
import scipy

## Somewhat arbitrary example
#np.random.seed(0)
#A = np.random.rand(4,4)

## Textbook example
#A = np.array([[2., 1., -1.], [-3., -1., 2.], [-2., 1., 2.]])

## A FP-sensitive example?
#A = scipy.linalg.hilbert(8)

# A pivoting example
A = np.array([[2., 1., -1.], [4., 2., -2.], [-2., 1., 2.]])

print("A = ")
print(A)

assert(A.shape[0] == A.shape[1])
n = A.shape[0]

E = np.empty((n,n))

# To compute the echelon form of A, we march across the columns of A.
# For diagonal point A[i,i]:
#   For row j >= 1:
#       E[i:,j] = A[i:,j] - (a[i,:] / a[i,i]) * A[i:,i]
#        
#   This sets column A[i:,i+1:] to zero, and modifies the rest.



# Compute the echelon (upper triangular) form of A

# This is like the "first term" in each element of E
E[:,:] = A[:,:]

for i in range(n-1):
    #print("i=", i)
    ## Collect the i column coefficient, weighted by the i diagonal
    #c = E[i+1:,i:i+1] / E[i,i]
    #print("c=")
    #print(c)

    ## Compute the diff row
    #r = E[i:i+1,i:]
    #print("r=")
    #print(r)
    
    ## Apply to all i+1 rows, which should force E[i+1:,i] to zero
    #E[i+1:,i:] = E[i+1:,i:] - c * r
    #print("E=")
    #print(E)

    ## Probably faster to do as a single line
    #E[i+1:,i:] = E[i+1:,i:] - (E[i+1:,i:i+1] / E[i,i]) * E[i:i+1,i:]

    # This expression also loops over the E[:i,:].
    #   These should all be zero, so it should be numerically equivalent.
    #   It may be faster, despite the extra work, due to pipelining and
    #   vectorization.
    E[i+1:,:] = E[i+1:,:] - (E[i+1:,i:i+1] / E[i,i]) * E[i:i+1,:]

# E is in echelon form
print("E = ")
print(E)


# Check the answer
def gaussian_elim(A):
    A = A.copy().astype(float)
    m, n = A.shape
    for i in range(min(m, n)):
        if A[i, i] == 0:
            continue  # skip zero pivot
        A[i+1:, :] -= (A[i+1:, i:i+1] / A[i, i]) * A[i:i+1, :]
    return A

E_check = gaussian_elim(A)

if np.allclose(E, E_check):
    print("Matches gaussin_elim()")
else:
    print("Does not match gaussian_elim()")
    print(E_check)
