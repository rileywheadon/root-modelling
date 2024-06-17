import json

with open("config.json") as f:
    CONFIG = json.load(f)
    PARAMS = CONFIG["PARAMS"]
    MODEL = CONFIG["MODEL"]
    
    rc = PARAMS["cells"][4]
    b = PARAMS["cells"][6]
    e = PARAMS["cells"][7]
    v = PARAMS["cells"][8]

    
def horizontal_connections(i, Ui0, Ui1):

    P, B = [], []
    column_i0 = get_column(i)
    column_i1 = get_column(i+1)
    j0, h0, j1, h1 = 0, 0, 0, 0
    
    while j0 < len(Ui0) and j1 < len(Ui1):
        
        b = min(h0 + Ui0[j0], h1 + Ui1[j1]) - max(h0, h1)
        
        # Add the background PIN to the list
        B.append((i, j0, i+1, j1, b))
        
        # Add right-facing PIN to the list (if the column type is correct)
        if (column_i0 in ["ER"] or j0 >= len(Ui0) - rc) and b > 0:
            P.append((i, j0, i+1, j1, b, "R"))
            
        # Add left-facing PIN to the list (if the column type is correct)
        if (column_i1 in ["EL"] or j0 >= len(Ui0) - rc) and b > 0:
            P.append((i+1, j1, i, j0, b, "L"))
            
        # Find the index of the higher cell, and increment it
        if h0 + Ui0[j0] < h1 + Ui1[j1]:
            h0 += Ui0[j0]
            j0 += 1
        else:
            h1 += Ui1[j1]
            j1 += 1

    return P, B


def vertical_connections(i, Ui):
    
    P, B = [], []
    column_type = get_column(i)
    
    for (j, u) in enumerate(Ui):
        
        # Check that the cell is not on the top edge
        # Then check for upward facing PINs
        
        if j != 0:
            B.append((i, j, i, j-1, PARAMS["width"]))
            
            if column_type in ["B"] or j >= len(Ui) - rc:
                P.append((i, j, i, j-1, PARAMS["width"], "U"))
            
        # Check that the cell is not on the bottom edge
        # Then check for downward facing PINs
        
        if j != len(Ui) - 1:
            B.append((i, j, i, j+1, PARAMS["width"]))
            
            if column_type in ["EL", "ER", "V"] or j >= len(Ui) - rc:
                P.append((i, j, i, j+1, PARAMS["width"], "D"))
    
    return P, B
     
    
def get_column(i):
    
    if i < b:
        return "B"
    elif i < b + e:
        return "ER"
    elif i < b + e + v:
        return "V"
    elif i < b + (2 * e) + v:
        return "EL"
    else:
        return "B"


def generate_pins(U):
    
    P, B = [], []
    
    # Define the list of PIN connections
    for i in range(len(U)):
        Pv, Bv = vertical_connections(i, U[i])
        P += Pv
        B += Bv
        
        if i != len(U) - 1:
            Ph, Bh = horizontal_connections(i, U[i], U[i+1])
            P += Ph
            B += Bh
            
    return P, B