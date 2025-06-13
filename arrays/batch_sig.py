from array_api_compat import get_namespace, to_device
import matplotlib.pyplot as plt
import numpy
import torch

def plot_batch_normalized_signal(batch):
    ns = get_namespace(batch)

    nbatch = batch.shape[0]

    norm = (batch - batch.mean(axis=1, keepdims=True)) / batch.std(axis=1, keepdims=True)

    norm_np = to_device(norm, 'cpu')
    
    plt.plot(norm_np.T)
    plt.show()


batch_numpy = numpy.random.randn(5, 100)
plot_batch_normalized_signal(batch_numpy)

batch_torch = torch.randn(5, 100)
plot_batch_normalized_signal(batch_torch)
