import numpy as np
import scipy

## Somewhat arbitrary example
#np.random.seed(0)
#A = np.random.rand(4,4)

# Textbook example
A = np.array([[2., 1., -1.], [-3., -1., 2.], [-2., 1., 2.]])

## A FP-sensitive example?
#A = scipy.linalg.hilbert(8)

## A pivoting example
#A = np.array([[2., 1., -1.], [4., 2., -2.], [-2., 1., 2.]])
A = np.array([[2., 1., -1.], [-2., 1., 2.], [4., 2., -2.]])

print("A = ")
print(A)

assert(A.shape[0] == A.shape[1])
n = A.shape[0]

E = np.empty((n,n))

# Copy A to E
#E[0,:] = A[0,:]

# To compute the echelon form of A, we march across the columns of A.
# For diagonal point A[i,i]:
#   For row j >= 1:
#       E[i:,j] = A[i:,j] - (a[i,:] / a[i,i]) * A[i:,i]
#        
#   This sets column A[i:,i+1:] to zero, and modifies the rest.

# Implementation notes:
# - The previous iteration sets E[:i,j] = 0, so we 
# - It may be numerically preferable to explicitly set E[i:,j] to zero?
# -     Does this help or hurt accuracy?

# Compute the echelon (upper triangular) form of A
E[:,:] = A[:,:]
#E[0,:] = A[0,:]
for i in range(1,n):
    c = E[i:,i-1:i] / E[i-1,i-1]
    E[i:,:] = E[i:,:] - c * E[i-1:i,:]

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
    print(E_check)
    print("Does not match gaussian_elim()")

