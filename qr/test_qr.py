import numpy as np
from qr import qr_decomp

# TODO: Fix the seed?

def test_decomp():
    A = np.random.rand(5,5)
    Q, R = qr_decomp(A)
    assert np.allclose(Q @ R, A)

def test_ortho():
    A = np.random.rand(5,5)
    Q, R = qr_decomp(A)
    I = np.eye(A.shape[0])
    assert np.allclose(Q.T @ Q, I)

def test_triu():
    A = np.random.rand(5,5)
    _, R = qr_decomp(A)
    assert np.allclose(R, np.triu(R))

def test_comp_numpy():
    A = np.random.rand(5,5)
    Q, R = qr_decomp(A)
    Q_np, R_np = np.linalg.qr(A)
    assert np.allclose(Q @ R, Q_np @ R_np)

# Not sure this will pass..
def test_rank_deficient():
    A = np.array([[1, 1], [2, 2], [3, 3]], dtype=float)
    Q, R = qr_decomp(A)
    assert np.allclose(Q @ R, A, atol=1e-12)

def test_identity():
    A = np.eye(5)
    Q, R = qr_decomp(A, positive_diagonals=True)
    I = np.eye(5)
    assert np.allclose(Q, I)
    assert np.allclose(R, I)
