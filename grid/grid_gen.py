import numpy as np
from math import sqrt

def semi_circle(shape, theta1, theta2, inner_radius, outer_radius):
    theta = np.linspace(theta1, theta2, shape[0])
    radius = np.linspace(inner_radius, outer_radius, shape[1])

    theta, radius = np.meshgrid(theta, radius, indexing='ij')

    return radius*np.cos(theta), radius*np.sin(theta)

def sine_channel(shape, start, end):
    x = np.linspace(start, end, shape[0])
    y = np.linspace(start, end, shape[1])

    x, y = np.meshgrid(x, y, indexing='ij')

    y += np.sin(x)

    return x, y

def bump(shape, domain_dims, radius):
    x, y = np.meshgrid(np.linspace(0, domain_dims[0], shape[0]), np.linspace(0, domain_dims[1], shape[1]), indexing='ij')

    x_loc = domain_dims[0]/2.

    for i in xrange(shape[0]):
        if abs(x[i, 0] - x_loc) < radius:
            ry = (radius**2 - (x[i, 0] - x_loc)**2)**0.5

            y[i, :] = np.linspace(ry, y[i, -1], shape[1])

    return x, y
