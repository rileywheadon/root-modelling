# Update Scheme

Simulation of the root proceeds in three steps:
1. A system of linear ODEs is created to represent the movement of each hormone between nodes. Then, a forward Euler method approximation is used to determine the new hormone concentrations. 
2. Next, the cells "grow" via the downward shifting of vertices in the model root. Hormone concentrations are adjusted to account for dilution caused by changes in cell volume. 
3. Finally, cells "divide" through the creation of new edges and vertices. Since the models assume hormones are distributed uniformly throughout the cell, daughter cells have an identical hormone profile to their parents.


