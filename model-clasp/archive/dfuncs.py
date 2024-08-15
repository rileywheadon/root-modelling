import json
import numpy as np
with open("data.json") as f:
    DATA = json.load(f)
    cpd, rot3, bes1, auxin, length, position = DATA.values()

def create_vector_function(params):
    xp, yp, degree, xlim, ylim = params.values()
    f = np.poly1d(np.polyfit(xp, yp, degree))
    return np.vectorize(lambda p : f(p) if p < xlim[1] else f(xlim[1]))

def dfuncs():
    return {
        "position_to_cpd": create_vector_function(cpd),
        "position_to_rot3": create_vector_function(rot3),
        "position_to_bes1": create_vector_function(bes1),
        "position_to_auxin": create_vector_function(auxin),
        "position_to_length": create_vector_function(length),
        "time_to_position": create_vector_function(position)
    }