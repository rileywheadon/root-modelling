# Transport Step

## Diffusion

All hormones are transported via diffusion between the nodes of the model.  The geometry of the model implies that every connection between nodes must pass through the cell wall, so the diffusion coefficient $D$ will be based on the permeability of the apoplast. The diffusion term for an arbitrary hormone $H$ in the $j$-th cell is:
$$D \sum_{i = 1}^{ N } b_{ij}(H_{i} - H_{j})$$

In this expression, $H_{j}$ represents the hormone concentration in the $j$-th cell, while $\{ H_{1}, \dots, H_{N} \}$ is the set of hormone concentrations in each of the cell's neighbours. The set $\{ b_{1j}, \dots, b_{Nj} \}$ contains the length of the boundary between the $j$-th node and every neighbour. For edge-vertex connections, this is assumed to be $\lambda$, the prescribed width of the cell wall.

## Biosynthesis and Degradation

Additionally, every hormone is subject to biosynthesis and degradation, which may vary throughout the root. Our model represents these processes in the term $s_{j} - d_{j}H_{j}$, where $s_{j}$ and $d_{j}$ represent the cell-specific biosynthesis and degradation rate respectively. 

## PIN-Mediated Transport

Auxin is also subject to unidirectional PIN-mediated efflux between the cytoplasm and cell wall. Our model makes the assumption that the PINs are distributed evenly along individual cell *edges*.  With this in mind, the PIN-mediated efflux term for a cell node $A_{j}$ and a edge node $A_{i}$ is $\pm A_{j}b_{ij}p_{ij}$ where $p_{ij}$ denotes the amount of PIN protein on the border. The term is positive in the edge since it is receiving auxin from the cell, and negative in the cell, since it is transferring auxin to the wall.

Additionally, every cell is assumed to have a small amount of background PIN protein that is expressed uniformly around the cell wall. This background efflux is set to $1/20^{ \text{th} }$ of the polar PIN efflux ([[@grieneisen2007]]).

## AUX1/LAX Influx Carriers

While PIN proteins determine the flux of auxin from the cell to cell wall, AUX1/LAX influx carriers define the flux from cell wall to cell. These influx carriers are distributed in a cell-type and zone specific manner, but are distributed uniformly around each cell. For a cell node $A_{j}$ and edge node $A_{i}$, influx due to AUX1/LAX is defined as  $\pm A_{i}b_{ij}I_{j}$ where $I_{j}$ denotes amount of AUX1/LAX in the cell.  This term is positive for the cell, since it is receiving auxin from the wall, and negative for the wall, since it is transferring auxin to the cell.

## Hormone Crosstalk

Although it has been proposed that brassinosteroid (BR) crosstalk with auxin varies spatially ([[@ackerman-lavert2020]]), our model is built off earlier research showing that BR is inhibited by auxin ([[@chaiwanon2015]]). To implement this crosstalk, we assume that BR degrades faster in regions with high auxin. Our model also implements the increased transcriptional regulation of PIN2 by BR ([[@hacham2012]]), and the inhibition of PIN1 and PIN3 by cytokinin ([[@ioio2008]]).

## Transport Equations

Within the cell, we use the following two differential equations for auxin transport/synthesis and brassinosteroid synthesis. Here, $A_{i}$ denotes the auxin within the cell, while $A_{j}$ denotes the auxin in an arbitrary wall region neighbouring the cell.
$$
\begin{align}
A_{i}' &= S_{A} - d_{A}A_{i} + \sum_{j} b_{ij} (IA_{j} - p_{ij}A_{i} - p_{bg}A_{i}) \\[5pt]
B' &= S_{B} - \left( d_{B} + d_{BA}\left( \frac{ A_{i}^{2} }{ K_{dBA}^{2} + A_{i}^{2} } \right)  \right) B
\end{align}
$$

Regions of the cell wall (edges) do not produce brassinosteroid, so $B' = 0$ for these nodes. However, edges still transmit auxin via PINs, so we use the following differential equation for auxin fluxes along on edges: 
$$
A_{j}' = \sum_{i} b_{ij}(P_{ij}A_{i} + P_{bg}A_{i} - IA_{j})
$$

PIN concentrations are also determined based on differential equations, where $C_{ki}$ and $B_{i}$ denote the cytokinin and brassinosteroid concentrations in the cell.
$$
P_{ij} = \begin{cases}
S_{p} - \left( d_{p} + d_{pC} \left( \dfrac{ C_{ki}^{2} }{ K_{dpC}^{2} + C_{ki}^{2} } \right)  \right) P_{ij} &\text{(PIN1, PIN3)} \\[5pt]
S_{p} - \left( d_{p}  - d_{pB}\left( \dfrac{ B_{i}^{2} }{ K_{dpB}^{2} + B_{i}^{2} } \right) \right)  P_{ij} &\text{(PIN2)}
\end{cases}
$$