# PyPFC
A python implementation of Phase Field Crystal solved with semi-implicit pseudospectral method


For the PFC, the conserved order parameter is a relative density $n(\boldsymbol{r},t)$ which evolves following a conservation equation as 
$$\frac{\partial n}{\partial t} = \nabla \cdot \left[ M \nabla \left( \frac{\delta F[n]}{\delta n} \right) \right]$$
with M being a mobility, $\delta F[n]/\delta n$ is the functional derivative of free-energy. The free-energy functional [^1] given by 
$$F[n(\boldsymbol{r},t)] = \int_V  \left[ \frac{1}{2} n \left(1+\nabla^2\right)^2 n + \frac{1}{4} n^2(2r+n^2) \right]\ d{\boldsymbol{r}}$$
with $r$ being a constant proportional to the temperature deviation from the melting point [^2].

## Example

The following figure is a result for the system with $M=1.0$, $r=-0.25$, $L_x=L_y=16\pi$, $N_x=N_y=2^8$, and $dt=0.1$. The initial condition is given by a normal distribution described by $n(x,y,t=0) = -0.285(1.0 - 0.02\mathcal{N}(0,1))$. The system evolves until $T=1500$.

### Plot of the last frame
![LastFrame](https://github.com/elvissoares/PyPFC/blob/master/pfc2d-crystal.png)

### Plot of the free-energy time evolution
![Freeenergy](https://github.com/elvissoares/PyPFC/blob/master/pfc2d-freenergy-crystal.png)

### An animated GIF
![GIF](https://github.com/elvissoares/PyPFC/blob/master/pfc2d-crystal.gif)

# References
[^1]: [Elder, K. R., and Martin Grant. Physical Review E 70, no. 5 (November 19, 2004): 051605](https://doi.org/10.1103/PhysRevE.70.051605)
[^2]:[Elder, K. R., Nikolas Provatas, Joel Berry, Peter Stefanovic, and Martin Grant. Physical Review B 75, no. 6 (February 14, 2007): 064107](https://doi.org/10.1103/PhysRevB.75.064107)
