# external libraries
import sys
sys.setrecursionlimit(1000000)
import os
import logging
import argparse

import math
import numpy as np
import torch
import matplotlib.pyplot as plt

sys.path.append("../optimizers/")

from optimizers.GD import *


def rosenbrock(tensor):
    x, y = tensor
    obj_val = (1 - x) ** 2 + 100 * (y - x ** 2) ** 2
    return obj_val


def plot_rosenbrok(iters):
    x = np.linspace(-2, 2, 250)
    y = np.linspace(-1, 3, 250)
    minimum = (1.0, 1.0)

    X, Y = np.meshgrid(x, y)
    Z = rosenbrock([X, Y])

    iter_x, iter_y = iters[0, :], iters[1, :]

    fig = plt.figure(figsize=(8, 8))

    ax = fig.add_subplot(1, 1, 1)
    ax.contour(X, Y, Z, 90, cmap='jet')
    ax.plot(iter_x, iter_y, color='r', marker='x')

    ax.set_title('Rosenbrock')
    plt.plot(*minimum, 'gD')
    plt.plot(iter_x[-1], iter_y[-1], 'rD')
    plt.savefig('docs/rosenbrock.png')



if __name__ == '__main__':

    path = "doc"
    try:
        os.mkdir(path)
    except OSError:
        print ("Directory %s already exists." % path)
    else:
        print ("Successfully created the directory %s " % path)


    initial_guess = (1.0, 0.5)
    tensor_size=2

    # training 
    x = torch.Tensor(initial_guess).requires_grad_(True)
    optimizer = GD([x], **{'lr': 1e-3, 'num_iter': 10})
    # with picture 
    # optimizer.optimize_history(rosenbrock, x, tensor_size, plot_rosenbrok)

    optimizer.optimize(rosenbrock, x)



