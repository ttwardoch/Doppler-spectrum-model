import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy
import random

from settings import *


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def unit_vector_at_angle(angle):
    """ Returns the unit vector at an angle from [1, 0] in degrees. clockwise  """
    angle = angle % 360
    return [np.cos(angle * np.pi / 180), np.sin(angle * np.pi / 180)]


def receiver_angle(angle):
    if 360 < angle <= 0:
        raise ValueError
    angle /= 2
    return np.cos(angle * np.pi / 180)


def cos_of_angle_between(v1, v2):
    """
    Returns the cos of angle between vectors 'v1' and 'v2'::
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)


def angle_between(v1, v2):
    """
    Returns the angle in radians between vectors 'v1' and 'v2'::
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def plot_3D(u, x, y, title=None):
    """Plot the latent variable field u given the list of x,y coordinates"""
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_trisurf(x, y, u, cmap='viridis', linewidth=0, antialiased=False)
    if title:  plt.title(title)
    plt.show()


def plot_2D(shifts, xi, yi, n, title=None, colors='viridis', crange=None):
    """Visualise count data given the index lists"""
    Z = -np.ones((n, n))
    for i in range(len(shifts)):
        Z[((i - i % n) // n, i % n)] = shifts[i]
    my_cmap = copy.copy(cm.get_cmap(colors))
    my_cmap.set_under('k', alpha=0)
    fig, ax = plt.subplots()
    if crange is None:
        crange = [-np.max(shifts), np.max(shifts)]
    im = ax.imshow(Z, origin='lower', cmap=my_cmap, clim=crange, extent=[0, PATCH_DIAMETER/1000, 0, PATCH_DIAMETER/1000])
    fig.colorbar(im)
    if title:  plt.title(title)
    plt.show(block=False)

def spread_from_bins(bins_labels, bins):
    dot = [] # Leave only labels that have anything in bins
    for label, bin in zip(bins_labels, bins):
        if bin > 0:
            dot.append(label)
    return max(dot) - min(dot)

def random_angles(random_type, ANGULAR_WIDTH_OF_PATCH):
    if random_type == "uniform":
        angle_perp = (random.random() - 0.5) * ANGULAR_WIDTH_OF_PATCH
        angle_tang = (random.random() - 0.5) * ANGULAR_WIDTH_OF_PATCH
    elif random_type == "radial":
        r = random.random() * ANGULAR_WIDTH_OF_PATCH / 2
        alpha = random.random() * np.pi * 2
        angle_perp = r * np.cos(alpha)
        angle_tang = r * np.sin(alpha)
    elif random_type == "gaussian":
        angle_perp = random.gauss(mu=0, sigma=STD_DEVIATION/2*ANGULAR_WIDTH_OF_PATCH)
        angle_tang = random.gauss(mu=0, sigma=STD_DEVIATION/2*ANGULAR_WIDTH_OF_PATCH)
    else:
        raise ValueError
    return angle_perp, angle_tang


def random_unit_circle_complex_value():
    theta = random.random()*2*np.pi
    x = np.cos(theta)
    y = np.sin(theta)
    return x + 1j*y
