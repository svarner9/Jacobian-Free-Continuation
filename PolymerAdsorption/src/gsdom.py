#============================================
#  Flory Solution near a Solid Surface
#     Example of Jacobian-free Pseudo-arclength
#     continuation for polymer adsorption
#
#   Written by Sam Varner and Chris Balzer (2024)
#   
#   See top-level README.md for citation.
#   Distributed under the MIT License.
#============================================

import os
import json
import numpy as np
from copy import deepcopy
from scipy.linalg import lstsq

class GroundStateDominance(object):
    # Constructor 
    def __init__(self, inputFile):
        # Default input parameters
        self.eta    = -0.05
        self.chi    = 2.0/3.0
        self.phiB0  = 1e-8
        self.dphiB  = 0.01
        self.N      = 100.0
        self.L      = 50.0
        self.dz     = 0.1
        self.steps  = 350
        self.mixType = 'anderson'
        self.anderM  = 8
        self.mixCoeff = 5e-3
        self.errorTol = 1.0e-6
        self.maxIters = 1e7
        self.continuation = True
        self.transform = 'sqrt'
        self.printDensity = False
        self.outputDir = './'

        # Read input file
        self.inputFile = 'input/' + inputFile
        self._read_input()

        # Calculate number of grid points
        self.M   = int(self.L/self.dz + 1)

        # Initialize density
        self.initialPhi()
        self.varphi = np.sqrt(self.phi)

        # Stepping
        self.stepNum = 0
        self.phiB = self.phiB0
        
        # Mixing
        self.iter = 0
        self.err = 1.0
       
        # Arrays for solvers
        self.numKeep = 10
        if self.continuation:
            self.numUnknowns = self.M + 1
        else:
            self.numUnknowns = self.M
        self.X        = np.zeros((self.numUnknowns,))
        self.F        = np.zeros((self.numUnknowns,))        
        self.X_Store  = np.zeros((self.numUnknowns,self.numKeep))
        self.F_Store  = np.zeros((self.numUnknowns,self.numKeep))

        # Anderson variables
        self.X_And  = np.zeros((self.numUnknowns,self.anderM))
        self.G_And  = np.zeros((self.numUnknowns,self.anderM))

        # Continuation variables
        self.s           = 0.0
        self.dS          = self.dphiB
        self.minDs       = 0.03
        self.maxDs       = 0.1
        self.order       = 0.0
        self.turnCount   = 0
        self.cont_vec    = np.zeros((self.steps))
        self.orderParam  = np.zeros((self.steps))
        self.s_vec       = np.zeros((self.numKeep))
        self.dYdS        = np.zeros((self.numUnknowns,1))
        self.dYdS_Store  = np.zeros((self.numUnknowns,self.numKeep))
        self._set_continuation()

        # Initialize X
        self._fill_X()

        # Data IO
        self.outputDir = 'output/' + self.outputDir + '/'
        self._set_output_file_structure()

    #===============================
    #   I/O Routines
    #===============================
    def _read_input(self):
        print(f"Reading input file: {self.inputFile} \n")
        with open(self.inputFile, 'r') as f:
            d = json.load(f)

        for key,value in d.items():
            if isinstance(value,str):
                value = value.lower()
            setattr(self,key,value)
        return d

    def _set_output_file_structure(self):
        if not os.path.exists(self.outputDir):
            os.makedirs(self.outputDir)

    def _writeParams(self):
        if self.stepNum == 0:
            tag = 'w'
        else:
            tag = 'a'
        with open(self.outputDir + 'params.dat',tag) as f:
            f.write(str(self.stepNum) + ' ' + f'{self.phiB:.5e}' + ' ' + f'{self._calc_orderParam():.5e}')
            if self.continuation:
                f.write(' ' + f'{self.s:.5f}' + ' ' + f'{self.turnCount:.0f}')
            f.write(' \n')

    def _fileWrite(method):
        def wrapper(self):
            fileName = method.__name__.split('_')[-1]
            file = self.outputDir + '/' + fileName  +  '_' + str(self.stepNum) + '.dat'
            method(self,file)
        return wrapper

    @_fileWrite
    def _write_density(self,fileName):
        if self.printDensity:
            with open(fileName,'w') as f:
                for i in range(self.M):
                    f.write(str(round(i*self.dz,3)) + ' ' + f'{self.phi[i]:.8e}'  + ' \n')


    #===============================
    #   Main functions
    #===============================
    def initialPhi(self):
        self.phi = self.phiB0*np.ones((self.M,))

    def solve(self):
        while (self.iter < self.maxIters) and (self.err > self.errorTol):
            self._update_field()
            self.iter += 1
        self._writeParams()
        self._write_density()

    def _update_field(self):
        self._fill_F()
        self._calc_error()
        self._anderson()
        self._get_X()

    def _fill_X(self):
        if self.transform == 'log':
            self.X[0:self.M] = np.log(self.phi)
        elif self.transform == 'sqrt':
            self.X[0:self.M] = np.sqrt(self.phi)
        else:
            self.X[0:self.M] = self.phi
        if self.continuation:
            self._set_continuation()

    def _get_X(self):
        if self.transform == 'log':
            self.phi = np.exp(self.X[0:self.M])
        elif self.transform == 'sqrt':
            self.phi = self.X[0:self.M]**2
        else:
            self.phi = self.X[0:self.M]
        self.varphi = np.sqrt(self.phi)
        if self.continuation:
            self._get_continuation()

    def _fill_F(self):
        chemPotential = 1.0/self.N * np.log(self.phiB/self.N) - np.log(1.0 - self.phiB) + self.chi*(1.0 - 2.0*self.phiB)
        for i in range(self.M):
            localChemPotential = 1.0/self.N * np.log(self.phi[i]/self.N) - np.log(1.0 - self.phi[i]) + self.chi*(1.0 - 2.0*self.phi[i])
            if i == 0:
                self.F[i] = -self.eta + 1.0/(6.0*self.varphi[i])*(-0.25*self.varphi[i+4] + 4.0/3.0*self.varphi[i+3] - 3.0*self.varphi[i+2] + 4.0*self.varphi[i+1] - 25.0/12.0*self.varphi[i])/self.dz
            elif i == self.M-1:
                self.F[i] = chemPotential - localChemPotential
            elif i == 1 or i == self.M-2:
                self.F[i] = chemPotential - localChemPotential + 1.0/(6.0*self.varphi[i])*(self.varphi[i+1] - 2.0*self.varphi[i] + self.varphi[i-1])/self.dz**2
            else:
                self.F[i] = chemPotential - localChemPotential + 1.0/(6.0*self.varphi[i])*(-self.varphi[i+2] + 16.0*self.varphi[i+1] - 30.0*self.varphi[i] + 16.0*self.varphi[i-1] - self.varphi[i-2])/12.0/self.dz**2
            
            if self.transform == 'log':
                self.F[i] *= 1.0
            elif self.transform == 'sqrt':
                self.F[i] *= 2.0*self.varphi[i]
        
        if self.continuation:
            if self.stepNum <= self.numKeep:
                self.F[-1] = 0.0
            elif self.continuation:
                self.F[-1] = self.dS - np.dot(self.dYdS,self.X-self.X0)

    def _calc_error(self):
        self.err = np.max(np.abs(self.F))

    def _picard(self):
        self.X += self.mixCoeff * self.F

    def _anderson(self):
        # Current column
        col = self.iter % self.anderM

        # Fill Anderson matrices
        self.X_And[:,col] = self.X
        self.G_And[:,col] = self.F

        if self.iter > self.anderM:
            try:
                # Temporary matrices
                tempX  = np.delete(self.X_And, col, axis=1) - self.X_And[:,col][:,None]
                tempG  = np.delete(self.G_And, col, axis=1) - self.G_And[:,col][:,None]
                lambdas,_,_,_ = lstsq(tempG,-self.G_And[:,col])

                self.X = (1 - self.mixCoeff)*(self.X + np.matmul(tempX,lambdas)) + self.mixCoeff*(self.X + self.F + np.matmul(tempX + tempG, lambdas) )
            except:
                self._get_X()
                self._write_density()
                exit()
        else:
            self._picard()

    def step(self):
        # Store values from step
        self.iter = 0
        self.err = 1
        self.stepNum += 1

        # Take step in parameter
        if not self.continuation:
            self.phiB = np.exp( np.log(self.phiB) +  self.dphiB )
        else:
            self._step_continuation()

    def _calc_orderParam(self): 
        temp1 = 0.0
        for i in range(0,self.M):
            temp = self.phi[i] - self.phiB
            if i == 0 or i == self.M-1:
                temp1 += 0.5 * temp
            else:
                temp1 += temp
        return temp1*self.dz
    
    def _get_continuation(self):
        self.phiB = np.exp(self.X[-1])
    
    def _set_continuation(self):
        self.X[-1] = np.log(self.phiB)

    def _step_continuation(self):
        # Convenient indices
        currentIndx = self.stepNum % self.numKeep
        prevIndx    = (self.stepNum - 1) % self.numKeep
        prevIndx2   = (self.stepNum - 2) % self.numKeep
        prevIndx3   = (self.stepNum - 3) % self.numKeep

        # Store previous profile
        self.s_vec[prevIndx]             = self.s
        self.cont_vec[self.stepNum-1]    = self.X[-1]
        self.orderParam[self.stepNum-1]  = self._calc_orderParam()
        self.X_Store[:,prevIndx]         = self.X
        self.X0                          = deepcopy(self.X)

        # Take predictor step
        if self.stepNum <= self.numKeep:
            # Natural continuation step
            self.X[-1]  += 0.03
            self.s      += 0.01
        else:
            # Parameters for gradient (non-uniform spacing)
            ds1 = self.s_vec[prevIndx]  - self.s_vec[prevIndx2]
            ds2 = self.s_vec[prevIndx2] - self.s_vec[prevIndx3]
            
            if self.stepNum >= 2*self.numKeep:
                c = (2*ds1 + ds2)/ds1/(ds1 + ds2)
                b = -(ds1 + ds2)/ds1/ds2
                a = ds1/ds2/(ds1 + ds2)
            else:
                c = 1.0/ds1
                b = -1.0/ds1
                a = 0.0

            # Calculate slope and norm
            self.dYdS = c*self.X_Store[:,prevIndx] + b*self.X_Store[:,prevIndx2]  + a*self.X_Store[:,prevIndx3]
            self.dYdS /= np.linalg.norm(self.dYdS)

            # Pick new ds (could use gradient to determine this)
            if np.abs(self.dYdS[-1]) > 0.8:
                self.dS *= 1.1
            else:
                self.dS *= 0.1

            if self.dS < self.minDs:
                self.dS = self.minDs
            elif self.dS > self.maxDs:
                self.dS = self.maxDs

            # Take step
            self.X  += self.dYdS*self.dS

            # Increment s
            self.s  += self.dS

           # Determine turn count
            if (self.X_Store[-1,prevIndx] - self.X_Store[-1,prevIndx2])*(self.X_Store[-1,prevIndx2] - self.X_Store[-1,prevIndx3]) < 0:
                self.turnCount += 1
                
        # Store continuation parameter
        self._get_continuation()
