import json

from computation.growth import grow
from computation.division import divide
from computation.auxin import flow

with open("config.json") as f:
    CONFIG = json.load(f)
    PARAMS = CONFIG["PARAMS"]
    MODEL = CONFIG["MODEL"]

    
def update_cells(network):
    
    # Growth step
    if MODEL["br"] == "direct" and not MODEL["static"]:
        network = grow(network, False)
        
    if MODEL["br"] == "clasp" and not MODEL["static"]:
        network = grow(network, True)

    # Division step
    if MODEL["auxin"] == "direct":
        network = divide(network, False)
        
    if MODEL["auxin"] == "pins":
        network = divide(network, True)
        network = flow(network)
    
    return network