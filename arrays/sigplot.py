from array_api_compat import get_namespace, to_device
import matplotlib.pyplot as plt
import numpy
import torch

a_np = numpy.sin(numpy.linspace(0, 10, 100))
a_torch = torch.sin(torch.linspace(0, 10, 100))

def plot_normalized_signal(a):
    ns = get_namespace(a)

    norm = (a - ns.mean(a)) / ns.std(a)
    
    # Although norm is a generic array, its contents may not be visible to the
    # CPU doing the plotting, so we need to move the data back to CPU.
    # This function tends to produce NumPy arrays.
    norm_np = to_device(norm, 'cpu')

    # Some trivia:
    # - to_device() will not copy norm if going from numpy to numpy.
    # - I think it's an Array API principle that positional arguments are
    #   rejected.  that's why the more common to_device(norm, device='cpu')
    #   doesn't work!

    # Anyway...now it can be safely plotted
    plt.plot(norm_np)
    plt.xlabel("Index")
    plt.ylabel("Normalized signal")
    plt.grid()
    plt.show()

plot_normalized_signal(a_np)
plot_normalized_signal(a_torch)
