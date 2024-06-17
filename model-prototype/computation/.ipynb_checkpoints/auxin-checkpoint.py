import json

from computation.network import get_column, generate_pins
from computation.boundary import boundary_influx, boundary_efflux

with open("config.json") as f:
    CONFIG = json.load(f)
    PARAMS = CONFIG["PARAMS"]
    MODEL = CONFIG["MODEL"]

    
# Moves auxin between cells via pins
def get_flow(U, C, P, B):
     
    # Create an empty dC (change in concentration) array with the same shape as C
    dC = []
    for Ci in C:
        dCi = []
        for c in Ci:
            dCi.append(0)
        dC.append(dCi)
        
    step = MODEL["step"]
    width = PARAMS["width"]
    strength = PARAMS["pin_strength"]
        
    # Iterate through P, adjusting dC as necessary
    for Pi in P:
        
        # Shifted auxin is (concentration * border * step)
        x0, y0, x1, y1, b, d = Pi
        shift = C[x0][y0] * b * step
        
        # Take away from (x0, y0), give to (x1, y1)
        dC[x0][y0] -= shift
        dC[x1][y1] += shift
        
    # Iterate through B, adjusting dC as necessary
    for Bi in B:
        
        x0, y0, x1, y1, b = Bi
        shift = ((C[x0][y0] - C[x1][y1]) * b * step) / strength
        
        dC[x0][y0] -= shift
        dC[x1][y1] += shift
    
    return dC


# Updates the PIN network and moves auxin between cells
def flow(network):
    
    U, C0, P0, B0, Z = network
    
    # Update the PIN network (if necessary)
    P1, B1 = generate_pins(U)

    # Move auxin between cells
    dC = get_flow(U, C0, P1, B1)
    dC = boundary_efflux(U, C0, dC)
    
    C1 = []
    for (C0i, dCi) in zip(C0, dC):
        C1i = [sum(c) for c in zip(C0i, dCi)]
        C1.append(C1i)

    # Set the boundary condition
    C1 = boundary_influx(C1)
    
    return (U, C1, P1, B1, Z)