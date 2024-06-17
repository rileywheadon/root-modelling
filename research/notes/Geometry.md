# Geometry

Our model represents the root through a graph of *nodes* and *connections*, where each node contains data about the hormones within it  ([[@band2014]]). Every node is either a *vertex*, an *edge*, or a *cell*. Cells are polygonal representations of the cytosol of a single root cell and the hormones within it. Surrounding each cell are edges of a constant width $W$, which represent regions of the cell wall. It is important to note that edges are also nodes which contain and transport hormones, albeit under a different set of rules than cells. Vertices represent the small, approximately square regions of the cell wall that connect edges to each other. For ease of computation, we take the area of each vertex to be $W^{2}$.

![[model-geometry.png|center|400]]