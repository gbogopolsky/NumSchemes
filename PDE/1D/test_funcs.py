import numpy as np

def gaussian(x, x0, sigma_x):
    """ Gaussian test function """
    return np.exp(- (x - x0)**2 / 2 / sigma_x**2)

def step(x, x0):
    """ Step test function """
    return np.where(abs(x - x0) < 0.5, 1.0, 0.0)

def packet_wave(x, x0, lam):
    """ Packet wave of spatial period lam """
    return np.where(abs(x - x0) < 0.5, np.sin(2 * np.pi / lam * (x - x0)), 0)
