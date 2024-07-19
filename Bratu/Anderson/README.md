## Bratu Problem (PAC with Anderson Acceleration)
This is an implementation of our PAC algorithm with Anderson acceleration where we solve the familiar Bratu problem with a second-order ODE that has instability. We use our PAC algorithm to trace out the full solution curve including thermal runaway.

The source code is located in ```src/bratu.py```. To run the code from within the current directory (```/Bratu/Anderson```), simply run the following in the command line:
```
$ python main.py --input input.json
```

### Input Parameters
There is one (optional) command line arguments that can be passed when running the code:
- ```--input```: The input file to use, always located in the ```input/``` folder (default: input.json)

Within the input file (```input/input.json```), there are a few parameters that can be changed:
```
{
    "steps": 500,
    "M": 50,
    "anderM": 20,
    "mixCoeff": 2e-5,
    "errorTol": 1e-6,
    "maxIters": 1e5
}
```
The meaning of the parameters is as follows:
- ```steps```: The maximum number of steps to take with PAC.
- ```M```: The number of gridpoints in each dimension (MxM grid).
- ```mixCoeff```: The mixing coefficient for Anderson Acceleration. Large values favor speed while small values favor stability. We recommend using values less than 5e-5.
- ```errorTol```: The error tolerance for the iterative solution at each corrector step.
- ```maxIters```: The maximum number of iterations for the iterative solution at each corrector step. Best practice is that each corrector step should reach `errorTol`, so this variable should be set to a large value in general.

### Output
There are two default outputs of the code. The first is the parameters file located at ```output/params.dat```. This file contains the information about the system at each step of the continuation process. The columns are:

0. Step number
1. Lambda parameter ($\lambda$)
2. Maximum temperature order parameter ($u_{max}$)
3. Number of iterations it took for corrector step.
4. Error of corrector step (should be less than ```errorTol```)
5. The pseudo-arclenth value ($s$)
6. The turn count. This is the number of times the continuation has turned around in the $\max(u)-\lambda$ plane. This is used to determine if there has been any turning points, and is helpful for plotting the different branches (stable vs unstable) in the continuation path.
7. Execution time for given PAC step (seconds).

The second is the iterations file located at ```output/iteration.dat```. This file contains the iteration information for the corrector step, and it is overwritten for each step of the PAC. It is mainly useful for checking to see if the calculation is stuck/not converging properly.

You can optionally output the profiles $u(x,y)$ for each step by adding the following argument to the input file:
```
"printDensity": "True"
```
