import json

import numpy as np

with open("config.json") as f:
    CONFIG = json.load(f)
    PARAMS = CONFIG["PARAMS"]
    MODEL = CONFIG["MODEL"]

    
# Takes a list of cells and returns a list of midpoints
def get_midpoints(Ui):
    ends = np.cumsum(Ui)
    starts = np.append(0, ends)[:-1]
    return (starts + ends) / 2


# Takes a column of cells and gets the differentiation zone size and gradient midpoint
def get_column_information(Ui):
    diff = 0
    length = 0
    cap = 0
    for j, u in enumerate(Ui):
        if u > PARAMS["differentiation"]:
            diff += u
        if j < len(Ui) - PARAMS["cells"][4]:
            cap += u
        length += u
    return (diff, length, cap)
    
    
# Takes a matrix of cells and gets the average differentiation zone size and gradient midpoint
def get_cell_information(U):
    
    rm, hm = PARAMS["cells"][2:4]
    
    diffs, lengths, caps = 0, 0, 0
    for Ui in U:
        d, l, c = get_column_information(Ui)
        diffs = max(diffs, d)
        lengths += l
        caps += c
        
    avg_cap = caps / len(U)
    avg_midpoint = (avg_cap + diffs) / 2
    return (diffs, avg_midpoint, avg_cap)
