import matplotlib.patches as patch

# Import the configuration file and set constants
import json
with open("config.json") as f:
    CONFIG = json.load(f)
    type_colors = CONFIG["DISPLAY"]["type_colors"]


class Hormones:
    def __init__(self):
        self.pins = {
            "1": 0,
            "2": 0,
            "3": 0,
            "4": 0,
            "7": 0
        }
        self.p = {
            "top": 0,
            "bottom": 0,
            "in": 0,
            "out": 0
        }
    
    def express_pins():
        return

class Vertex:
    def __init__(self, x: int, y: int):
        self.x = x
        self.y = y

    def __repr__(self):
        return f"({self.x}, {self.y})"

    def __eq__(self, other):
        if isinstance(other, Vertex):
            return (self.x == other.x and self.y == other.y)
        return False

class Wall:
    def __init__(self, v1: Vertex, v2: Vertex):
        
        # Vertex a must have lower x, (or lower y to break ties)
        verts = sorted([v1, v2], key = lambda v : (v.x, v.y))
        self.v1, self.v2 = verts[0], verts[1]

        # Initialize the neighbours list and hormones
        self.neighbours = []
        self.hormones = Hormones()

    def __repr__(self):
        return f"{self.v1.__repr__()}<->{self.v2.__repr__()}" 

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, Wall):
            return (self.v1 == other.v1 and self.v2 == other.v2)
        return False
    
    def add_neighbour(self, n):
        if n not in self.neighbours:
            self.neighbours.append(n)
            n.add_neighbour(self)
    
    def view(self):
        return [self.v1.x, self.v2.x], [self.v1.y, self.v2.y]

class Cell:
    def __init__(self, V: [Vertex], t: str):
        self.V = V
        self.type = t
        self.neighbours = []
        self.hormones = Hormones()

    def add_neighbour(self, n):
        if n not in self.neighbours:
            self.neighbours.append(n)
            n.add_neighbour(self)
            
    def height(self):
        ys = [v.y for v in self.V]
        return max(ys) - min(ys)

    def view(self):
        xy = [(v.x, v.y) for v in self.V]
        fill = type_colors[self.type]
        return patch.Polygon(xy, ec="gray", lw = 0.5, fc=fill)