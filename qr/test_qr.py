import numpy as np
from qr import qr_decomp
#import torch
#from qr_array import qr_decomp

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

def test_overdetermined():
    A = np.random.rand(8,4)
    Q, R = qr_decomp(A)
    # maybe i need param testing?
    I = np.eye(A.shape[0])
    assert np.allclose(Q.T @ Q, I)
    assert np.allclose(R, np.triu(R))
    assert np.allclose(Q @ R, A)

def test_underdetermined():
    A = np.random.rand(4,8)
    Q, R = qr_decomp(A)
    # maybe i need param testing?
    I = np.eye(A.shape[0])
    assert np.allclose(Q.T @ Q, I)
    assert np.allclose(R, np.triu(R))
    assert np.allclose(Q @ R, A)

def test_reduced():
    A = np.random.rand(4,8)
    Q, R = qr_decomp(A, mode='reduced')
    assert np.allclose(Q @ R, A)

def test_positive_diagonals():
    A = np.random.rand(5,5)
    _, R = qr_decomp(A, positive_diagonals=True)
    assert np.all(np.diag(R) > 0.)

def test_identity():
    A = np.eye(5)
    Q, R = qr_decomp(A, positive_diagonals=True)
    I = np.eye(5)
    assert np.allclose(Q, I)
    assert np.allclose(R, I)

def test_zero_column():
    A = np.random.rand(5,5)
    A[:,2] = 0.
    Q, R = qr_decomp(A)
    assert np.allclose(Q @ R, A)

#def test_torch():
#    A = torch.rand(5,5)
#    Q, R = qr_decomp(A)
#    torch.allclose(Q @ R, A)
