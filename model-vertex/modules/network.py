import numpy as np
from modules.classes import *

# Import the configuration file and set constants
import json
with open("config.json") as f:
    CONFIG = json.load(f)

    # Column parameters
    columns = CONFIG["COLUMNS"]["types"]
    column_widths = CONFIG["COLUMNS"]["widths"]
    max_size = CONFIG["COLUMNS"]["max_height"]
    mz_size = CONFIG["COLUMNS"]["mz_height"]
    min_size = CONFIG["COLUMNS"]["min_height"]
    ez_rows = CONFIG["COLUMNS"]["ez_rows"]
    mz_rows = CONFIG["COLUMNS"]["mz_rows"]

    # Cap parameters
    cap_rows = CONFIG["CAP"]["types"]
    cap_columns = CONFIG["CAP"]["columns"]
    min_cap_size = CONFIG["CAP"]["min_height"]
    max_cap_size = CONFIG["CAP"]["max_height"]

    # Calculated parameters
    rows = ez_rows + mz_rows
    width = sum([column_widths[t] for t in columns])
    cap_width = width + (cap_rows.count("OUT") * column_widths["OUT"] * 2)


# Generate cell columns
def get_height(i):
    if i < ez_rows:
        step = (max_size - mz_size) / ez_rows
        return max_size - (i * step)
    return mz_size 

def generate_columns():
    ws = [0] + [column_widths[t] for t in columns]
    hs = [0] + [get_height(i) for i in range(rows)]
    xs = np.cumsum(ws) - width / 2
    ys = np.cumsum(hs)
    
    C = []
    for x0, x1, t in zip(xs, xs[1:], columns):
        for y0, y1 in zip(ys, ys[1:]): 
            V = [Vertex(x0, y0), Vertex(x1, y0), Vertex(x1, y1), Vertex(x0, y1)]
            C.append(Cell(V, t))
            
    return C, max(ys)


# Generate the root cap
def interpolate(s, h, base):

    # Define the function for the curve of vertices
    f = lambda x : (h / s**3) * (s**3 - np.abs(x ** 3)) + base
    xs = np.linspace(-s, s, 100)
    ys = f(xs)

    # Get the total length of the curve ts
    dx, dy = xs[1:]-xs[:-1], ys[1:]-ys[:-1]
    ds = np.array((0, *np.sqrt(dx*dx + dy*dy)))
    ts = np.cumsum(ds)
    
    # Interpolate
    xinter = np.interp(np.linspace(0, ts[-1], cap_columns + 1), ts, xs)
    yinter = np.interp(np.linspace(0, ts[-1], cap_columns + 1), ts, ys)
    return list(zip(xinter, yinter))

def generate_curves(base):
    ws = [0] + [column_widths[t] for t in cap_rows]
    ss = (cap_width / 2) - np.cumsum(ws)
    hs = np.cumsum(np.linspace(min_cap_size, max_cap_size, 6))[::-1]
    return [interpolate(s, h, base) for s, h in zip(ss, hs)]


def generate_cap(base):
    C, curves = [], generate_curves(base)
    for row0, row1 in zip(curves, curves[1:]):
        for j, verts in enumerate(zip(row0, row0[1:], row1[1:], row1)):
            t = "LAT" if (j < 3 or j > 6) else "COL"
            V = [Vertex(v[0], v[1]) for v in verts]
            C.append(Cell(V, t))
    return C, curves[-1]
    

# Generate the QC
def generate_qc(base, rim):
    C, s = [], column_widths["VAS"]
    T = [[0, 1, 2, 3], [3, 4, 5], [5, 6, 7], [7, 8, 9, 10]]
    B = [[-2*s, -3*s], [0, -s, -2*s], [2*s, s, 0], [3*s, 2*s]]
    for Ti, Bi in zip(T, B):
        v0 = [Vertex(rim[t][0], rim[t][1]) for t in Ti]
        v1 = [Vertex(b, base) for b in Bi]
        C.append(Cell(v0 + v1, "QUI"))
    return C



# Gets the edges for a cell, subdividing if necessary
def get_edges(cell):
    
    # Determine the number of subdivisions. Cap cells have no divisions.
    divisions = int(cell.height() // 16)
    if (cell.type in ["LAT", "COL", "QUI"]) or divisions == 0:
        divisions = 1
        
    # Get the top and bottom edges
    verts = sorted(cell.V, key = lambda v : v.y)
    middle = len(verts) // 2
    t = sorted(verts[:middle], key = lambda v : v.x)
    b = sorted(verts[middle:], key = lambda v : v.x)
    
    # Get the left and right edges points for the new cell(s)
    lx = [t[0].x - ((t[0].x - b[0].x) * (i / divisions)) for i in range(divisions + 1)]
    ly = [t[0].y - ((t[0].y - b[0].y) * (i / divisions)) for i in range(divisions + 1)]
    left = [Vertex(x, y) for x, y in zip(lx, ly)][1:-1]
    rx = [t[-1].x - ((t[-1].x - b[-1].x) * (i / divisions)) for i in range(divisions + 1)]
    ry = [t[-1].y - ((t[-1].y - b[-1].y) * (i / divisions)) for i in range(divisions + 1)]
    right = [Vertex(x, y) for x, y in zip(rx, ry)][1:-1]
        
    # Generate the wall segments for the new cells
    w = []
    outer = t + right + b[::-1] + left[::-1]
    for v0, v1 in zip(outer, outer[1:] + [outer[0]]):
        w.append(Wall(v0, v1))
    
    # Create the layers of the new cells
    layers = [t] + [[l, r] for l, r in zip(left, right)] + [b]
        
    # Iterate through the layers and create the new cells
    c = []
    for l0, l1 in zip(layers, layers[1:]):
        c.append(Cell(l0 + l1[::-1], cell.type))
    
    return c, w

# Generates a list of edges from a list of cells
def generate_edges(cells):
    
    # Generate the list of cell walls
    C, W = [], []
    for c in cells:
        c, w = get_edges(c)
        C += c
        W += w

    # Remove any duplicates
    W = list(set(W))
    
    # Add wall/cell connections to the network
    for w in W:
        for c in C:
            
            # Check each cell to see if contains both vertices of w
            if w.v1 in c.V and w.v2 in c.V:
                c.add_neighbour(w)
                
    return W, C


# Generates the entire network of cells and edges
def generate_network():
    columns, base = generate_columns()
    cap, rim = generate_cap(base)
    qc = generate_qc(base, rim)
    C = columns + cap + qc
    return generate_edges(C)