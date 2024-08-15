import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from computation.gradient import gradient, get_phi

with open("config.json") as f:
    CONFIG = json.load(f)
    PARAMS = CONFIG["PARAMS"]
    MODEL = CONFIG["MODEL"]
    VISUAL = CONFIG["VISUAL"]

    
# Returns patches for a column of cells
def draw_column(x, dx, Ui, Ci):
    
    # Get the number of rows in the cap
    rc = PARAMS["cells"][4]
    
    # Initialize the three collections (differentiated, regular, cap)
    y, pid, cid, pir, cir, pic, cic = 0, [], [], [], [], [], []
    
    # Add each cell to the correct collection
    for j, (u, c) in enumerate(zip(Ui, Ci)):
        
        box = plt.Rectangle((x, y), dx, u)
        if j >= len(Ui) - rc:
            pic.append(box)
            cic.append(c)
        elif u > PARAMS["differentiation"]:
            pid.append(box)
            cid.append(c)
        else:
            pir.append(box)
            cir.append(c)
            
        y += u
    
    return (pid, cid, pir, cir, pic, cic)


# Draws a matrix of cells, returns the patch collection
def draw_cells(U, C, dx):

    # Plot the cells
    x, pd, cd, pr, cr, pc, cc = 0, [], [], [], [], [], []
    
    for Ui, Ci in zip(U, C):
        pid, cid, pir, cir, pic, cic = draw_column(x, dx, Ui, Ci)
        pd += pid
        cd += cid
        pr += pir
        cr += cir
        pc += pic
        cc += cic
        x += dx
    
    # Add the patches
    cells0 = PatchCollection(pd, edgecolor="black", cmap=VISUAL["cmap"], hatch="x")
    cells1 = PatchCollection(pr, edgecolor="black", cmap=VISUAL["cmap"])
    cells2 = PatchCollection(pc, edgecolor="darkgray", cmap=VISUAL["cmap"])
    
    # Set colourmap arrays and limits
    cells0.set_array(cd)
    cells1.set_array(cr)
    cells2.set_array(cc)
    
    cells0.set_clim(0, VISUAL["auxin_max"])
    cells1.set_clim(0, VISUAL["auxin_max"])
    cells2.set_clim(0, VISUAL["auxin_max"])
    
    return (cells0, cells1, cells2)


# Adds an individual arrow to the patch collection
def arrow(U, Pi, xs):
    
    # Get the aspect ratio and unpack the PIN information
    aspect = VISUAL["height"] / 3
    x0, y0, x1, y1, b, d = Pi
    
    # Find the heights of the bottom / top of both cells to align arrows
    h0, h1 = sum(U[x0][:y0]), sum(U[x1][:y1])
    u0, u1 = U[x0][y0], U[x1][y1]
    
    # Set the arrow size and position for horizontal and vertical arrows
    if d in ["L", "R"]:
        py = (max(h0, h1) + min(h0 + u0, h1 + u1)) / 2
        dx = (x1 - x0) * (xs / 3)
        dy = 0
    if d in ["U", "D"]:
        px = (x0 * xs) + (xs / 2)
        dx = 0
        dy = (y1 - y0) * (xs / 3) * aspect
    
    # Adjust for arrow direction
    if d == "L":
        px = (x0 * xs) + (1 * xs) / 6
    if d == "R":
        px = (x0 * xs) + (5 * xs) / 6
    if d == "U":
        py = h0 + (xs * aspect) / 6
    if d == "D":
        py = (h0 + u0) - (xs * aspect) / 6
        
    # Adjust arrow dimensions based on the aspect ratio
    w = 0.01 if dx == 0 else 0.01 * aspect
    hw = 0.04 if dx == 0 else 0.04 * aspect
    hl = 0.02 * aspect if dx == 0 else 0.02
    
    return patches.FancyArrow(px, py, dx, dy, width = w, head_width = hw, head_length = hl)
    

# Creates a collection of all arrows
def draw_pins(U, P, xs):  
    patches = [arrow(U, pin, xs) for pin in P]
    return PatchCollection(patches, facecolor = VISUAL["pin_colour"])


# Draws a the colorbar on the main axes
def draw_colorbar(ax, data):
    
    # Draw the colorbar
    ins = inset_axes(
        ax,
        width = "5%",
        height = "30%",
        loc = "upper right",
        borderpad = 4,
    )
    
    ins.set(title = "Auxin")
    fig = ax.figure
    fig.colorbar(data, cax = ins, ticks = [0, 0.2, 0.4, 0.6, 0.8, 1])
    
    return ax
    

# Draws zone lines on the main axes
def draw_zones(ax, Z):
    
    # Unpack Z
    e, m, c, rc = Z
    
    # Draw the zone lines
    zone1 = ax.axhline(y = e, color = 'b', linestyle = '-', label = "DZ / EZ") 
    zone2 = ax.axhline(y = m, color = 'm', linestyle = '-', label = "EZ / MZ")
    zone3 = ax.axhline(y = c, color = 'r', linestyle = '-', label = "MZ / Cap")
    ax.legend(loc = "upper left")
    
    return ax


# Draws the main axes
def axes_main(ax, U, C, P, Z):
    
    # Determine cell width constants
    dx = 1 / len(U)
    
    # Get cells
    cells0, cells1, cells2 = draw_cells(U, C, dx)
    
    # Add the zone borders
    ax.clear()
    if VISUAL["show_zones"]:
        ax = draw_zones(ax, Z)
    
    # Configure the plot
    ax.add_collection(cells0)
    ax.add_collection(cells1)
    ax.add_collection(cells2)
    ax.set_box_aspect(1)
    ax.set(xlim=(-1, 2), ylim=(VISUAL["height"], 0))
    ax.set_xticks([])
    
    # If PINs are unused, return as is
    if MODEL["auxin"] != "pins":
        return ax
    
    # Otherwise, draw PINs and add colorbar
    if VISUAL["show_pins"]:
        pins = draw_pins(U, P, dx)
        ax.add_collection(pins)
        
    return draw_colorbar(ax, cells0)


# Draws the secondary gradient axes
def axes_gradient(ax, U, Z):
    
    # Generate the gradient, phi, and find the midpoint
    e, m, c, rc = Z
    br, auxin = gradient(PARAMS["gradient"], m)
    phi = get_phi(PARAMS["phi"])
    
    # Draw the axes, based on the settings in MODEL
    ax.clear()
    
    ys = np.linspace(1, VISUAL["height"], 200)
    
    if MODEL["br"] != "none":
        br_plot, = ax.plot(br(ys), ys, label = "BR", color = "orange")
    if MODEL["br"] == "clasp":
        phi_plot, = ax.plot(phi(br(ys)), ys, label = r"$\phi$", color = "green")
    if MODEL["auxin"] == "direct":
        auxin_plot, = ax.plot(auxin(ys), ys, label = "Auxin", color = "purple")
    
    ax.tick_params(axis = 'y', which = 'both', left = False)
    ax.legend(loc = "upper left")
    
    return ax
