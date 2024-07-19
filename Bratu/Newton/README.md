## Bratu Problem (Continuation with Newton Solver)
This is an implementation of Newton-based PAC solver for the familiar Bratu problem. With this code you can solve for the full solution curve including the unstable (runaway) branch.

The source code is located in ```src/newton.py```. To run the code from within the current directory (```/Bratu/Newton```), simply run the following in the command line:
```
$ python main.py -N 50 -dS 0.5
```
### Input Parameters
There are two (optional) command line arguments that can be passed when running the code:
- ```-N```: The number of gridpoints in each dimension (default: 50)
- ```-dS```: The stepsize for arc-length continuation (default: 0.5)

### Output
There is a single output of the code, which is the parameters file located at ```output/params.dat```. This file contains the information about the system at each step of the continuation process. The columns are:

0. Step number
1. Lambda parameter ($\lambda$)
2. Maximum temperature order parameter ($u_{max}$)
3. The pseudo-arclenth value ($s$)
4. The turn count. This is the number of times the continuation has turned around in the $\max(u)-\lambda$ plane. This is used to determine if there has been any turning points, and is helpful for plotting the different branches (stable vs unstable) in the continuation path.
