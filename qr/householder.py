import numpy as np

debug = False

# Householder implementation of Q R factorization

print("A = ")
A = np.random.rand(4,4)
print(A)

# R is not yet triangular, but it will be!
R = np.array(A)
Q = np.eye(*A.shape)

ncols = A.shape[1]
x = np.zeros(ncols)
e = np.zeros(ncols)
v = np.zeros(ncols)

# What's up with the final iteration?  Sometimes the matrix blows up
for c in range(ncols-1):
    # Copy the next column vector of the submatrix into x
    x[:c] = 0.
    x[c:] = R[c:,c]

    # Construct the projection basis vector
    e[:] = 0.
    e[c] = 1.

    # Flush the contents of the previous iteration
    v[:c] = 0.
    
    # Compute the hyperplane unit vector v:
    #   x + x_proj = V
    #   V = x + |x| e
    # NOTE: There are two sign conventions, depending on whether the hyperplane
    #   is based on V = x + x_proj or x = x_proj + V.  (NumPy uses the first.)
    v = x + np.linalg.norm(x) * e
    v = v / np.linalg.norm(v)

    # Construct the Householder reflection operation Hx = x - 2 <x.v> v
    H = np.eye(*R.shape) - 2. * np.outer(v, v)

    # Update Q and R
    #   R is now "more" triangular, and another H is applied to Q.
    # A = Q R = Q (H R') = (Q H') R'
    Q = Q @ H
    R = H @ R

    if debug or c == ncols-2:
        print("Q = ")
        print(Q)

        print("R = ")
        print(R)


# And check the result:
print("Is Q orthogonal?", np.allclose(Q.T @ Q, np.eye(*Q.shape), atol=1e-15))
print("Is R upper triangular?", np.all(np.tril(R,k=-1) < 1e-15))
print("Does A = QR?", np.allclose(Q @ R, A, atol=1e-15))
