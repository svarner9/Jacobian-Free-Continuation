## Bratu Equation (Liouville–Bratu–Gelfand)
We provide code to solve for the solution curve of the Bratu Equation. The equation can be used to model the temperature in a reactor where the walls are at a fixed temperature. For an exothermic reaction, there is a feedback loop where the reaction acts as a heat source to further promote more reaction, which then releases more heat. However, due to the fixed temperature of the boundary, this heat can be dissipated to keep the system under control. If the reaction is very exothermic, then the heat flux at the boundary is not sufficient to modulate the temperature, and the reaction will runaway. The simplest form of the equation is a second-order ODE:

$$\nabla^2 u + \lambda e^u = 0 $$

where $u$ is the temperature, and $\lambda>0$ is a measure of how exothermic the reaction is. In the paper, and in this code, we solve the equation in 2d on the domain $x\in[0,1], y\in[0,1]$. We use the dirichlet boundary conditions of $u=0$ everywhere on the boundary of the square domain.

As explained in the paper, there are both stable and unstable solutions to this equation. In order to capture the unstable solutions, we must start in the stable region, and take steps along the solution curve using an arclength continuation technique.

### Provided Codes/Methods

Here, we provide two different codes that can compute the full solution curve. The code located in `Anderson` contains our Anderson-based PAC algorithm (outlined in the paper). The code located in `Newton` contains a traditional Newton-based PAC algorithm that relies on the Jacobian.

To run the provided codes, change directories into the respective folders and following the instructions in the `README.md` file.

```
$ cd Anderson

or

$ cd Newton
```
