# ⭐️ Root Growth Project

>[!definition]+ Links
>- [[⭐️ Root Growth Notes]]
>- [[📔 Research Journal]]
>- [[❓ Open Questions]]

>[!theorem]+ Reading List
>- [ ] [[@marconi2021]]

>[!remark]+ Development Todo (Single Cell)
> 

>[!remark]+ Development Todo (2D)
>
>**GOAL**: Build model directly from ODEs
>- [ ] Define edge locations (top, bottom, in, out). `Cell.connections` should be a tuple `(direction, Wall)`
>- [ ] Write code to visualize the different hormone levels within the root (take a hormone as a parameter, return a plot) 
>- [ ] Think about how to define a model from ODEs
>	- `Hormones` object stores PIN1, PIN2, PIN3
>	- `Hormones` has `express_pins` which converts PIN1, PIN2, PIN3 to `p_top`, `p_in`, `p_out`, `p_bottom`
>--- 
>- [ ]  Prescribe the [[Initial Conditions]] for every variable:
>- [ ] Code the transport step of the update scheme
>	- Iterate through the list of edges, cells, nodes
>	- Create a list of "orders" (changes in each hormone for each node)
>	- Iterate through the list of orders and apply to each node
>- [ ] Code the growth step of the update scheme
>- [ ] Code the division step of the update scheme
>	- Also divide subcellular regions that are too big ($30 \mu m$)

>[!example]+ Writing Todo
>- [ ] Write [[Ethylene]]
>- [ ] Write [[PLETHORA Crosstalk Models]]
>- [ ] Write [[Cytokinin Crosstalk Models]]
>- [ ] Write [[Ethylene Crosstalk Models]]
