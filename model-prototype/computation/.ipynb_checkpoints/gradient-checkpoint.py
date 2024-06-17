import numpy as np


# Generates a gradient (tuple of logistic functions) with provided slope and midpoint
def gradient(slope, midpoint):
    br = lambda x : 1 / (1 + np.exp(slope * (x - midpoint)))
    auxin = lambda x : 1 - (1 / (1 + np.exp(slope * (x - midpoint))))
    return (br, auxin)


# Returns a function that takes the level of BR in [0, 1] and returns the MT alignment in [0, 1]
def get_phi(params):
    s, t, u, v, w, k_in0, k_inB, k_on, k_off, C_m, n = params
    
    # Compute intermediate variables a, b, S, T, U
    a = lambda x : (k_on * x) / (k_inB + k_off)
    b = lambda x : (k_in0 + (k_on * x) - (k_off * a(x))) / (v * w)
    S = lambda x : u * a(x) * b(x)
    T = lambda x : b(x) - ((u * a(x)) / w)
    U = lambda x : 0 - (1 / w) - (s / t)
    
    # Compute steady state values of R0, RB, C
    R0 = lambda x : ((0 - T(x)) + np.sqrt((T(x) ** 2) - (4 * S(x) * U(x)))) / (2 * S(x))
    RB = lambda x : (a(x) * R0(x))
    C = lambda x : (b(x) * R0(x)) - (1 / w)
    
    # Return the function phi
    return lambda x : (C_m ** n) / ((C_m ** n) + (C(x) ** n))