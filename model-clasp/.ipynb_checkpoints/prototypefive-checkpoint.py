import numpy as np
import pandas as pd
from numba import njit

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
    return BR / np.max(BR)


# Data parsing function for mutant data
def parse_mutant_data(prefix, raw, radius):

    # Filter by the columns that contain the specified prefix
    columns = [c for c in raw.columns if c.startswith(prefix)]
    df = (raw.loc[:, columns])
    
    # Determine the cell lengths using a cylindrical assumption
    df_lengths = (df
        .sub(2 * np.pi * radius ** 2)
        .div(2 * np.pi * radius)
    )

    # Compute the mean cell length for each cell number
    df_means = pd.DataFrame({c : df_lengths.mean(axis = 1) for c in columns})
    
    # Determine the cell positions, filling in the mean when necessary
    df_positions = (df_lengths
        .fillna(df_means)
        .cumsum(axis = 0)
    )
    
    # Transform data into (position, length) tuples
    df_tuples = (pd
        .concat([df_positions, df_lengths])
        .stack()
        .groupby(level=[0,1])
        .apply(np.array)
        .unstack()
        .to_numpy()
        .flatten()
        .tolist()
    )
    
    # Filter out NaN length values and negative length values
    mask = lambda x : len(x) == 2 and x[1] > 0
    data = np.array(list(filter(mask, df_tuples)))
    return data, len(columns)