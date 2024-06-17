# Model Planning

- Parameters (can vary spatially)
- Hormones (these move between cells)
- Intracellular factors (PINS, growth factors, etc.)
- Transport equations / crosstalk (System of DEs)

## Ideas

Only cells and edges, one node per cell.

Need a way to generate an initial distribution of spatially-varied parameter values. Perhaps a 2D array with the values in each region:

```python
[[epi/mat, cor/mat, end/mat, vas/mat],
 [epi/elo, cor/elo, end/elo, vas/elo],
 [epi/div, cor/div, end/div, vas/div],
 [lat/cap, col/cap]]
```

Hormones are stored in a `Hormones` object in each cell. Use `setattr()` to build this object from the config.
- `auxin`, `br`, `pin1`, `pin2`, `pin3`

Crosstalk is defined using `lambda` functions, which can apply to either cell-to-edge or edge-to-cell transport relationships.

Going to need an `express_pins()` function that takes `Edge, Cell` and determines the requested `p_ij` for that cell-edge connection.
