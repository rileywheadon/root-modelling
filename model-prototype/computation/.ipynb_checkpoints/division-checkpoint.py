import json

import numpy as np

from computation.gradient import gradient
from computation.zones import get_midpoints

with open("config.json") as f:
    CONFIG = json.load(f)
    PARAMS = CONFIG["PARAMS"]
    MODEL = CONFIG["MODEL"]


# Uses a random number generator to divide the cell
def divide_column(U0i, C0i, auxin, pins, m):
    
    M = get_midpoints(U0i)
    P = auxin(M)
    
    U1i, C1i, y = [], [], 0
    for j, (u, c, p) in enumerate(zip(U0i, C0i, P)):
        
        y += u
        
        # Only allow cell division in:
        # - Cells below the elongation zone
        # - Cells not in the root cap
        # - Cells smaller than the division threshold
        below_e = y > m
        not_cap = j < len(U0i) - PARAMS["cells"][4]
        not_small = u > PARAMS["division"]
        
        # Check randomly for cell division
        threshold = c if pins else p
        threshold *= MODEL["step"] * PARAMS["division_scale"]
        r = np.random.default_rng().random()
        division = r < threshold
        
        if below_e and not_cap and not_small and division:
            U1i.append(u / 2)
            U1i.append(u / 2)
            C1i.append(c)
            C1i.append(c)
        else:
            U1i.append(u)
            C1i.append(c)
    
    return U1i, C1i


# Divides cells based on a computed auxin gradient
def divide(network, pins):
    
    if MODEL["static"]:
        return network
    
    U0, C0, P, B, Z = network
    e, m, c, rc = Z
    br, auxin = gradient(PARAMS["gradient"], m)

    # Divide the cells
    U1, C1 = [], []
    for U0i, C0i in zip(U0, C0):
        U1i, C1i = divide_column(U0i, C0i, auxin, pins, m)
        U1.append(U1i)
        C1.append(C1i)
        
    return (U1, C1, P, B, Z)
