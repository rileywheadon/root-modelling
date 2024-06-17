# Model Technicals

### Computation

A single function `update` takes all the data and returns a new set of data. Within it, there should be sequential calls to
- `update_pins` - Determines the new PIN list
- `update_auxin` - PIN-mediated auxin flow
- `update_cells` - Cell division and growth

### Visualization

`initialize_plot` returns the starting figure and a list of axes. All animation is done through the `draw` functions below.

`axis_main` draws the main axes (the cells). It will call:
- `draw_zones` - Draws the zone boundaries
- `draw_cells` - Draws the cells
- `draw_pins`- Draws PINs on top of the cells

`axis_gradient` draws the BR / auxin / phi gradient axis

`axis_colorbar` draws the colorbar axis for auxin levels

### Animation

1. Make a single call to `initialize_plot`
2. Iterate over `update` and construct an array `animationInfo` that contains all of the data needed to produce the animation.
3. Declare an `update_frame` function that takes an `info` from `animationInfo` and returns the new axes
4. Call `animate` and export the resulting file

