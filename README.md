# PyPFC
A python implementation of Phase Field Crystal solved with semi-implicit pseudospectral method


For the PFC, the conserved order parameter is a relative density $n(\boldsymbol{r},t)$ evolves following a conservation equation as 
$$\frac{\partial n}{\partial t} = \nabla \cdot \left[ M \nabla \left( \frac{\delta F[n]}{\delta n} \right) \right]$$
with M being a mobility, $\delta F[n]/\delta n$ is the functional derivative of free-energy. The free-energy functional given by 
$$F[n(\boldsymbol{r},t)] = \int_V  \left[ \frac{1}{2} n \left(1+\nabla^2\right)^2 n + \frac{1}{4} n^2(2r+n^2) \right]\ d{\boldsymbol{r}}$$
with $r$ being a constant proportional to the temperature deviation from the melting point \cite{Elder2004,Elder2007}.

## Example

The following figure is a result for the system with $M=1.0$,$r=-0.25$, $L=16\pi$, $N=2^8$, and $dt=0.1$. The initial condition is given by a normal distribution described by $n(x,y,t=0) = -0.285(1.0 - 0.02\mathcal{N}(0,1))$.

### Plot of the last frame
![LastFrame](https://github.com/elvissoares/PyPFC/blob/master/pfc2d-crystal.png)

### Plot of the free-energy time evolution
![Freeenergy](https://github.com/elvissoares/PyPFC/blob/master/pfc2d-freenergy-crystal.png)

### An animated GIF
![GIF](https://github.com/elvissoares/PyPFC/blob/master/pfc2d-crystal.gif)
