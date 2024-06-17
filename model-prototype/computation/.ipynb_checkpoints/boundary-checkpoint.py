import json

import numpy as np

from computation.network import get_column

with open("config.json") as f:
    CONFIG = json.load(f)
    PARAMS = CONFIG["PARAMS"]
    MODEL = CONFIG["MODEL"]
    
    step = MODEL["step"]
    width = PARAMS["width"]
    strength = PARAMS["pin_strength"]
    influx = PARAMS["influx"]

    
def boundary_influx(C):
    
    for i in range(len(C)):
        C[i][0] = influx
        
    return C


def boundary_efflux(U, C, dC):

    # Left edge auxin loss
    for j, c in enumerate(C[0]):
        dC[0][j] -= (C[0][j] * U[0][j] * step) / strength
        
    # Right edge auxin loss
    xmax = len(C) - 1
    for j, c in enumerate(C[xmax]):
        dC[xmax][j] -= (C[xmax][j] * U[xmax][j] * step) / strength
        
    # Bottom edge auxin loss
    for i, Ci in enumerate(C):
        ymax = len(Ci) - 1
        dC[i][ymax] -= (C[i][ymax] * width * step) / strength

    # Handle outflow through the top border cells
    for i in range(len(C)):
        if get_column(i) == "B":
            dC[i][0] -= C[i][0] * width * step
    
    return dC
