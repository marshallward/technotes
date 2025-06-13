from array_api_compat import get_namespace
import numpy as np
import torch

def mean_and_std(a):
    ns = get_namespace(a)
    return ns.mean(a), ns.std(a)

a_np = np.array([1, 2, 3, 4, 5])
a_torch = torch.tensor([1, 2, 3, 4, 5], dtype=torch.float32)

print(mean_and_std(a_np))
print(mean_and_std(a_torch))
