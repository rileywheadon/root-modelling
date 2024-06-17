# Growth Step

The current model approaches the problem of symplastic growth (see [[Growing Root Models]]) using a similar approach to Salvi et al. (2020) ([[@salvi2020]]). In order to ensure all cells at a given height have the same growth rate, growth factors are averaged across lateral subsections of the root, and applied uniformly to all cells in that subsection. To accommodate a biologically accurate root cap, cells in the root cap are fixed. Additionally, cells stop growing once they have reached a certain size, as in the model proposed in Pavelescu et al. (2018) ([[@pavelescu2018]]).

## CLASP Model

Let $C$ be the CLASP, $R_{0}$ the empty receptors, an $R_{B}$ the bound receptors. $B$ denotes the quantity of brassinosteroid in the cell. The interactions between these variables observed *in vivo* can be modelled using the following system of ordinary differential equations (ODEs):

$$
\begin{cases}
C' = c_{\text{max}} - p(1 + qR_{B})C \\[5pt]
R_{0}' = u(1 + vC) - k_{\text{on}}R_{0}B- k_{\text{in}}^{ 0 }R_{0}  + k_{\text{off}}R_{B} \\[5pt]
R_{B}' = k_{\text{on}}R_{0}B - k_{\text{in}}^{ B }R_{B} - k_{\text{off}}R_{B}
\end{cases}
$$

CLASP production $C'$ attains a maximum $c_{max}$ when $C = 0$ within the cell. The parameters $p$ and $q$ modulate the restriction of CLASP production due to the current CLASP level $C$ and number of bound receptors $R_{B}$.

Empty receptors $R_{0}$ are produced at a rate $u(1 + vC)$ with respect to the CLASP $C$ and the parameters $u$ and $v$.  Empty receptors bind to  BR at a rate of $k_{\text{on}}R_{0}B$, while bound receptors unbind at a rate of $k_{\text{off}}R_{B}$. Additionally, both bound and unbound receptors are removed from the at rates of $k_{in}^{ B }R_{B}$ and $k_{in}^{ 0 }R_{0}$ respectively. 

## Time Scale Analysis

We assume that the intracellular regulation of CLASP takes place on a shorter time scale compared to the processes of cellular growth and division. Therefore, the levels of CLASP, empty receptors, and bound receptors maintains a steady state relative to BR in the context of the entire root. This allows us to make a quasi-steady state approximation of our system of ODEs, where $C' = R_{0}' = R_{B}' = 0$. Solving the resulting nonlinear system of equations yields a function $C(B)$ that defines CLASP in terms of BR.

As discussed previously, longitudinal root growth is requires organized transverse arrays of microtubules within the growing cells. When CLASP is high, microtubules instead form transfacial bundles, which restrict growth. Since CLASP is inversely proportional to growth, we define the growth rate of a cell $\phi \in [0, 1]$ to be given by the logistic function:
$$
\phi (C) = \frac{ C_{m}^{ n } }{ C_{m}^{ n } + C^{ n } }
$$

Here, $C_{m}$ governs the location of the transition between high and low growth, while $n$ governs how gradual or sudden the aforementioned transition is. By composing $\phi$ with the function $C(B)$ described above, we obtain $\phi(B)$, a single function that defines the growth of a cell with respect to its brassinosteroid concentration.