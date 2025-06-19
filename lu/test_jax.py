import jax.numpy as jnp
from array_api_compat import to_device
from lu_array import LU

A = jnp.array([[2.0, 3.0], [4.0, 5.0]])
A_api = A.__array_namespace__().array(A)  # or just use A

P, L, U = LU(A_api)
PA = jnp.matmul(P, A)
LU_ = jnp.matmul(L, U)

print("max error", jnp.max(jnp.abs(PA - LU_)))

