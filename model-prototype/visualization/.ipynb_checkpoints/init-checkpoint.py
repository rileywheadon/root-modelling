import json

import matplotlib.pyplot as plt

from visualization.axes import axes_main, axes_gradient
from computation.network import generate_pins
from computation.boundary import boundary_influx 

with open("config.json") as f:
    CONFIG = json.load(f)
    PARAMS = CONFIG["PARAMS"]
    MODEL = CONFIG["MODEL"]
    

# Generates the initial arrangement of cells and concentrations
def generate_network():
    
    re, he, rm, hm, rc, hc, b, e, v = PARAMS["cells"]
    rows = re + rm + rc
    cols = (2 * b) + (2 * e) + v
    U, C = [], []
    
    # Initialize the U, C arrays
    for i in range(cols):
        
        Ui, Ci = [], []
        for j in range(rows):
            if j < re: 
                Ui.append(he)
            elif j < re + rm: 
                Ui.append(hm)
            else: 
                Ui.append(hc)
            Ci.append(0)
            
        U.append(Ui)
        C.append(Ci)
        
    # Construct the Z tuple
    e = 0
    m = re * he
    c = m + (rm * hm)
    Z = (e, m, c, rc)
        
    # Return U, C, P if PINs are off
    if MODEL["auxin"] != "pins":
        return (U, C, [], [], Z)
    
    # Set the initial condition and PIN concentrations
    C = boundary_influx(C)
    P, B = generate_pins(U)

    return (U, C, P, B, Z)


# Draws a new plot
def initialize_plot(network):
    
    # Unpack the network tuple
    U, C, P, B, Z = network
    
    # Create the figure and axes
    fig, (ax1, ax2) = plt.subplots(
        nrows = 1, 
        ncols = 2, 
        sharey = 'row', 
        figsize = (18, 13),
        gridspec_kw = {'width_ratios': [3, 1]},
        constrained_layout = True
    )
    
    # Draw the main axis
    ax1 = axes_main(ax1, U, C, P, Z)

    # Draw the gradient
    ax2 = axes_gradient(ax2, U, Z)
    
    return (fig, ax1, ax2)