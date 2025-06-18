import numpy as np
import scipy

from lu_array import LU

## Somewhat arbitrary example
np.random.seed(0)
A = np.random.rand(4,4)

# Textbook example
B = np.array([[2., 1., -1.], [-3., -1., 2.], [-2., 1., 2.]])

# A high condition number example
C = scipy.linalg.hilbert(4)

# A pivoting example
D = np.array([[2., 1., -1.], [4., 2., -2.], [-2., 1., 2.]])

for M in A, B, C, D:
#for M in (A,):
    print("M = ")
    print(M)

    P, L, U = LU(M)

    # E is in echelon form
    print("U = ")
    print(U)
    print("L = ")
    print(L)
    print("P = ")
    print(P)
    
    
    # Check the answer
    
    if np.allclose(P @ M, L @ U):
        print("P M = L U")
    else:
        print("ERROR! P M /= L U!")
        print("M=")
        print(M)
        print("P @ L @ U =")
        print(P.T @ L @ U)

    print(40*"-")
