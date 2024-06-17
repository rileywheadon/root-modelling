import json

from computation.gradient import gradient, get_phi
from computation.zones import get_midpoints, get_cell_information

with open("config.json") as f:
    CONFIG = json.load(f)
    PARAMS = CONFIG["PARAMS"]
    MODEL = CONFIG["MODEL"]
    
    
def compose(f, g):
    return lambda *a, **kw: f(g(*a, **kw))


def grow(network, use_phi):
    
    U0, C, P, B, Z0 = network
    e0, m0, c0, rc = Z0
    br, auxin = gradient(PARAMS["gradient"], m0)
    phi = get_phi(PARAMS["phi"])
    vfunc = compose(phi, br) if use_phi else br
    
    # Grow the cells
    U1 = U0
    J = [0 for Ui in U0]
    Y = [Ui[0] for Ui in U0]
    y, step = min(Y), min(Y)
    
    while True:
        
        # Check that we are not in root cap or differentiation zone
        if y > e0 and y < c0 and step != 0:
        
            # Apply growth
            dy = vfunc(y) * step * MODEL["step"] * PARAMS["growth_scale"]
            for Ui, j in zip(U1, J):
                Ui[j] += dy

        # Increment arrays
        minY = min(Y)
        idx = Y.index(minY)
        J[idx] += 1

        # Exit condition
        if len(U0[idx]) == J[idx]: 
            break

        # Otherwise move to next subsection
        Y[idx] += U0[idx][J[idx]]
        step = minY - y
        y = minY
    
    # Update Z
    e1, m1, c1 = get_cell_information(U1)
    Z1 = (e1, m1, c1, rc)
    
    return (U1, C, P, B, Z1)