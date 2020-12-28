import numpy as np
from numba import njit

@njit(cache=True)
def advance(res, u, sigma, scheme):
    """ Temporal advancement of 1D finite difference numerical schemes """
    if scheme == 'FOU':
        res[0] = sigma * (u[0] - u[-1])
        for i in range(1, len(u)):
            res[i] = sigma * (u[i] - u[i - 1])
    elif scheme == 'LW':
        res[0] = sigma / 2 * (u[1] - u[-2]) - sigma**2 / 2 * (u[1] + u[-2] - 2 * u[0])
        res[-1] = sigma / 2 * (u[1] - u[-2]) - sigma**2 / 2 * (u[1] + u[-2] - 2 * u[-1])
        for i in range(1, len(u) - 1):
            res[i] = (sigma / 2 * (u[i + 1] - u[i - 1]) 
                        - sigma**2 / 2 * (u[i + 1] + u[i - 1] - 2 * u[i]))
    elif scheme == 'SOU':
        res[0] = (sigma / 2 * (3 * u[0] - 4 * u[-2] + u[-3]) - sigma**2 / 2 * (u[-3] - 2 * u[-2] + u[0]))
        res[1] = (sigma / 2 * (3 * u[1] - 4 * u[0] + u[-2]) - sigma**2 / 2 * (u[-2] - 2 * u[0] + u[1]))
        for i in range(2, len(u)):
            res[i] = (sigma / 2 * (3 * u[i] - 4 * u[i - 1] + u[i - 2]) 
                        - sigma**2 / 2 * (u[i - 2] - 2 * u[i - 1] + u[i]))
    elif scheme == 'TOU':
        res[0] = (sigma / 6 * (11 * u[0] - 18 * u[-1] + 9 * u[-2] - 2 * u[-3])
                    - sigma**2 / 2 * (2 * u[0] - 5 * u[-1] + 4 * u[-2] - u[-3]))
        res[1] = (sigma / 6 * (11 * u[1] - 18 * u[0] + 9 * u[-1] - 2 * u[-2])
                    - sigma**2 / 2 * (2 * u[1] - 5 * u[0] + 4 * u[-1] - u[-2]))
        res[2] = (sigma / 6 * (11 * u[2] - 18 * u[1] + 9 * u[0] - 2 * u[-1])
                    - sigma**2 / 2 * (2 * u[2] - 5 * u[1] + 4 * u[0] - u[-1]))
        for i in range(3, len(u)):
            res[i] = (sigma / 6 * (11 * u[i] - 18 * u[i - 1] + 9 * u[i - 2] - 2 * u[i - 3])
                    - sigma**2 / 2 * (2 * u[i] - 5 * u[i - 1] + 4 * u[i - 2] - u[i - 3]))
    elif scheme == 'C1':
        res[0] = - (0.5 * sigma * (sigma + 1) * u[-2]
                        + (1 - sigma**2) * u[0]
                        + 0.5 * sigma * (sigma - 1) * u[1] - u[0])
        res[-1] = - (0.5 * sigma * (sigma + 1) * u[-2]
                        + (1 - sigma**2) * u[-1]
                        + 0.5 * sigma * (sigma - 1) * u[1] - u[-1])
        for i in range(1, len(u) - 1):
            res[i] = - (0.5 * sigma * (sigma + 1) * u[i - 1]
                        + (1 - sigma**2) * u[i]
                        + 0.5 * sigma * (sigma - 1) * u[i + 1] - u[i])
    elif scheme == 'C2':
        res[0] = - ((sigma - 1) * sigma * (sigma + 1) * (sigma + 2) / 24 * u[-3]
                - (sigma - 2) * sigma * (sigma + 1) * (sigma + 2) / 6 * u[-2]
                + (sigma - 2) * (sigma - 1) * (sigma + 1) * (sigma + 2) / 4 * u[0]
                - (sigma - 2) * (sigma - 1) * sigma * (sigma + 2) / 6 * u[1]
                + (sigma - 2) * (sigma - 1) * sigma * (sigma + 1) / 24 * u[2]) + u[0]
        res[-1] = - ((sigma - 1) * sigma * (sigma + 1) * (sigma + 2) / 24 * u[-3]
                - (sigma - 2) * sigma * (sigma + 1) * (sigma + 2) / 6 * u[-2]
                + (sigma - 2) * (sigma - 1) * (sigma + 1) * (sigma + 2) / 4 * u[-1]
                - (sigma - 2) * (sigma - 1) * sigma * (sigma + 2) / 6 * u[1]
                + (sigma - 2) * (sigma - 1) * sigma * (sigma + 1) / 24 * u[2]) + u[-1]
        res[1] = - ((sigma - 1) * sigma * (sigma + 1) * (sigma + 2) / 24 * u[-2]
                - (sigma - 2) * sigma * (sigma + 1) * (sigma + 2) / 6 * u[-1]
                + (sigma - 2) * (sigma - 1) * (sigma + 1) * (sigma + 2) / 4 * u[1]
                - (sigma - 2) * (sigma - 1) * sigma * (sigma + 2) / 6 * u[2]
                + (sigma - 2) * (sigma - 1) * sigma * (sigma + 1) / 24 * u[3]) + u[1]
        res[-2] = - ((sigma - 1) * sigma * (sigma + 1) * (sigma + 2) / 24 * u[-4]
                - (sigma - 2) * sigma * (sigma + 1) * (sigma + 2) / 6 * u[-3]
                + (sigma - 2) * (sigma - 1) * (sigma + 1) * (sigma + 2) / 4 * u[-2]
                - (sigma - 2) * (sigma - 1) * sigma * (sigma + 2) / 6 * u[0]
                + (sigma - 2) * (sigma - 1) * sigma * (sigma + 1) / 24 * u[1]) + u[-2]
        for i in range(2, len(u) - 2):
            res[i] = - ((sigma - 1) * sigma * (sigma + 1) * (sigma + 2) / 24 * u[i - 2]
                    - (sigma - 2) * sigma * (sigma + 1) * (sigma + 2) / 6 * u[i - 1]
                    + (sigma - 2) * (sigma - 1) * (sigma + 1) * (sigma + 2) / 4 * u[i]
                    - (sigma - 2) * (sigma - 1) * sigma * (sigma + 2) / 6 * u[i + 1]
                    + (sigma - 2) * (sigma - 1) * sigma * (sigma + 1) / 24 * u[i + 2]) + u[i]
        
@njit(cache=True)
def iterations(nt, res, u, sigma, scheme):
    # iterations
    for _ in range(nt):
        advance(res, u, sigma, scheme)
        u -= res