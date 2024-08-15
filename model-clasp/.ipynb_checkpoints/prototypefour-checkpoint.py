from scipy.stats.qmc import Sobol
from numba import njit
import json
import numpy as np

# Compute the moving average of BR precursor data with a window of size n 
# Prescribe n = 50, since this is a low sensitivity parameter
@njit
def moving_average(data, n = 50): 
    
    output = np.zeros(1000)
    counts = np.zeros(1000)

    for row in data:
        left = max(0, int(row[0]) - n)
        right = min(1000, int(row[0]) + n)
        output[left:right] = output[left:right] + row[1]
        counts[left:right] = counts[left:right] + 1

    return np.divide(output, counts)


# Get the BR signalling (a weighted average of CPD and ROT3)
# Prescribe bias = 0.5, since this is a low sensitivity parameter
# Also transform by a factor of 100 to prevent numerical issues
@njit
def get_br(vP, CPD, ROT3, bias = 0.5):

    # Get the CPD and ROT3 data in terms of time by interpolating with vP
    ps = np.linspace(0, 999, 1000)
    CPD_MA = np.interp(vP, ps, moving_average(CPD))
    ROT3_MA = np.interp(vP, ps, moving_average(ROT3))
    BR = (CPD_MA * bias) + ROT3_MA * (1 - bias)
    return BR / 100

# Performs a Quasi-Monte Carlo Sensitivity Analysis
#   m = log2(N), where N is the number of sample points
#   d is the number of parameters (dimensionality)
#   bounds is an array of shape (d, 2)
#   cost is the cost function to use for analysis
def quasi_monte_carlo(m, d, bounds, cost):
    N = 2 ** m
    sbl = Sobol(d)
    
    # Generate random N x d matrices A, B
    A = sbl.random_base2(m)
    B = sbl.random_base2(m)

    # Build A_B matrices from A and B
    AB = []
    for i in range(d):
        ABi = np.copy(A)
        ABi[:, i] = B[:, i]
        AB.append(ABi)
        
    AB = np.array(AB)

    # Set the bounds, then transform the unit hypercube to the parameter space
    @njit
    def transform(row):
        # bounds = [(10, 200), (0, 1), (0, 10), (0, 0.1), (0, 0.1)]
        return [b[0] + (b[1] - b[0]) * p for p, b in zip(row, bounds)] 

    # Run the model on a single matrix M
    @njit
    def simulate_matrix(M):
        return np.array([cost(transform(row)) for row in M])

    # Run the model on all (d + 2) matrices
    fA = simulate_matrix(A)
    fB = simulate_matrix(B)
    fAB = np.array([simulate_matrix(ABi) for ABi in AB])

    # Get the variance of the output value Y
    varY = np.var(np.concatenate((fA, fB)))

    # Get the first order sensitivity index
    @njit
    def first_order_sensitivity(fA, fB, fABi):
        return (1 / N) * np.dot(fB, np.abs(fABi - fA)) / varY

    # Get the structural contribution of sensitivity
    @njit
    def structural_sensitivity(fB):
        return (1 / N) * np.sum(np.square(fB)) / varY

    # Get the total effect index
    @njit
    def total_effect_index(fA, fABi):
        return (1 / 2 * N) * np.sum(np.square(fA - fABi)) / varY

    # Compute the sensitivity and total effect indices
    sensitivities = [first_order_sensitivity(fA, fB, fABi) for fABi in fAB]
    structural = [structural_sensitivity(fB)]
    total_effects = [total_effect_index(fA, fABi) for fABi in fAB]

    # Return the sensitivity and total effect indices
    return sensitivities, total_effects