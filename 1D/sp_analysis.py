#!/Users/cheng/code/envs/dl/bin/python
import os
import numpy as np
import cmath
import matplotlib.pyplot as plt


def G_LW(phi, sigma):
    """ Return the Lax-Wendroff scheme amplification factor """
    return 1 + sigma**2 * (np.cos(phi) - 1) - 1j * sigma * np.sin(phi)


def G_FOU(phi, sigma):
    """ Return the First Order Upwind scheme amplification factor """
    return 1 - sigma * (1 - np.exp(-1j * phi))


def G_SOU(phi, sigma):
    """ Return the Second Order Upwind scheme amplification factor """
    return (1 - sigma / 2 * (3 - 4 * np.exp(-1j * phi) + np.exp(-2j * phi))
                + sigma**2 / 2 * (1 - 2 * np.exp(-1j * phi) + np.exp(-2j * phi)))

def G_TOU(phi, sigma):
    """ Return the Third Order Upwind scheme amplification factor -> unstable at low CFL """
    return (1 - sigma / 6 * (11 - 18 * np.exp(-1j * phi) + 9 * np.exp(-2j * phi) - 2 * np.exp(-3j * phi))
                    + sigma**2 / 2 * (2 - 5 * np.exp(-1j * phi) + 4 * np.exp(-2j * phi) - np.exp(-3j * phi)))

def G_C1(phi, sigma):
    return (0.5 * sigma * (sigma + 1) * np.exp(-1j * phi)
                        + (1 - sigma**2)
                        + 0.5 * sigma * (sigma - 1) * np.exp(1j * phi))

def G_C2(phi, sigma):
    """ Fourth order centered stable scheme """
    return ((sigma - 1) * sigma * (sigma + 1) * (sigma + 2) / 24 * np.exp(-2j * phi)
                    - (sigma - 2) * sigma * (sigma + 1) * (sigma + 2) / 6 * np.exp(-1j * phi)
                    + (sigma - 2) * (sigma - 1) * (sigma + 1) * (sigma + 2) / 4
                    - (sigma - 2) * (sigma - 1) * sigma * (sigma + 2) / 6 * np.exp(1j * phi)
                    + (sigma - 2) * (sigma - 1) * sigma * (sigma + 1) / 24 * np.exp(2j * phi))

def errors(G, phi, sigma):
    """ Computation of diffusion and dispersion error for a constant advection
    speed problem """
    G_num = G(phi, sigma)
    diff_err = abs(G_num)
    disp_err = np.zeros_like(phi)
    disp_err[0] = 1
    disp_err[1:] = np.array([- cmath.phase(G_num[i]) / sigma / phi[i] for i in range(1, len(phi))])
    return diff_err, disp_err
