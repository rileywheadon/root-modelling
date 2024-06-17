# Model Data Structures


`network = (U, C, P, Z)` is the data structure that contains all cells
	
`U` is the ragged 2D array of cell heights
`C` is the ragged 2D array of cell auxin concentrations

We use `U` to generate `P` (PINS), a list with the following format
```python
[(x0, y0, x1, y1, b, d), ...]
```

- `x0, y0` is the initial cell in `C` where auxin is coming from
- `x1, y1` is the cell in `C` where auxin is travelling to
- `b` denotes the size of the boundary between the two cells
- `d` is the direction of the PIN (left, right, up, down)

We also generate `B` (Background PINS), a list with this format:
```python
[(x0, y0, x1, y1, b)]
```

`Z = (e, m, c, rc)` is a tuple of zone information
- `e` is the start of the elongation zone
- `m` is the start of the meristematic zone
- `c` is the start of the root cap
- `rc` is the number of cells in the root cap

