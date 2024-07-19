#============================================
#  2D Bratu Problem
#       Solved using Jacobian-free PAC
#
#   Written by Sam Varner and Chris Balzer (2024)
#   
#   See top-level README.md for citation.
#   Distributed under the MIT License.
#============================================

import os
import json
import numpy as np
from time import time
from copy import deepcopy
from scipy.linalg import lstsq


class Bratu(object):
    def __init__(self, inputFile):

        # Default input file parameters
        self.steps = 500          # Maximum number of contination steps
        self.M = 50               # Number of grid points in each dimension
        self.anderM  = 30         # Number of histories for Anderson acc.
        self.mixCoeff = 1e-5      # Mixing coefficient
        self.errorTol = 1e-6      # Error tolerance for corrector step
        self.maxIters = int(1e5)  # Max number of iters for corrector step
        
        # Default continuation parameters
        self.continuation = True
        self.printDensity = False
        self.outputDir = './'

        # Read input file
        self.inputFile = 'input/' + inputFile
        self._read_input()

        # Create grid
        self.xVec = np.linspace(0,1,self.M)
        self.yVec = np.linspace(0,1,self.M)
        self.h = self.xVec[1] - self.xVec[0]
        self.M = self.M - 2

        # Initialize phi (same as u vector)
        self.initialPhi()
        
        # Stepping
        self.stepNum = 0
        self.lamb = 0.0
        self.lambVec = []
        self.orderParamVec = []
        
        # Mixing
        self.iter = 0
        self.anderIter = 0
        self.err = 1.0
        self.mix0 = self.mixCoeff
       
        # Arrays for solvers
        self.numKeep = 5
        if self.continuation:
            self.numUnknowns = self.M**2 + 1
        else:
            self.numUnknowns = self.M**2
        self.X        = np.zeros((self.numUnknowns,))
        self.F        = np.zeros((self.numUnknowns,))  
        self.F2d      = np.zeros((self.M,self.M))      
        self.X_Store  = np.zeros((self.numUnknowns,self.numKeep))
        self.F_Store  = np.zeros((self.numUnknowns,self.numKeep))

        # Anderson variables
        self.X_And  = np.zeros((self.numUnknowns,self.anderM))
        self.G_And  = np.zeros((self.numUnknowns,self.anderM))

        # Continuation variables
        self.s           = 0.0
        self.dS          = 0.1
        self.minDs       = 1.0
        self.maxDs       = 3.0
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
            f.write(str(self.stepNum) + ' ' + f'{self.lamb:.5e}' + ' ' + f'{self._calc_orderParam():.5e}' + ' ' + f'{self.iter:.0f}' + ' ' + f'{self.err:.5e}')
            if self.continuation:
                f.write(' ' + f'{self.s:.5f}' + ' ' + f'{self.turnCount:.0f}' + ' ' + f'{self.simTime:.5f}')
            f.write(' \n')

    def _fileWrite(method):
        def wrapper(self):
            fileName = method.__name__.split('_')[-1]
            file = self.outputDir + fileName  +  '_' + str(self.stepNum) + '.dat'
            method(self,file)
        return wrapper

    @_fileWrite
    def _write_density(self,fileName):
        if self.printDensity:
            with open(fileName,'w') as f:
                for i in range(self.M):
                    for j in range(self.M):
                        f.write(str(round(self.xVec[j+1],3)) + ' ' + str(round(self.yVec[i+1],3)) + ' ' + f'{self.phi[i][j]:.8e}' + ' \n')
                    # f.write(str(round((i+1)*self.h,3)) + ' ' + f'{self.phi[i]:.8e}'  + ' \n')

    # Make one file that get overwritten each step
    def _write_iteration(self,fileName='iteration.dat'):
        if self.iter == 0:
            tag = 'w'
        else:
            tag = 'a'
        with open(self.outputDir + fileName, tag) as f:
            f.write(str(self.iter) + ' ' + f'{self.err:.5e}' + '\n')

    def initialPhi(self):
        self.phi = np.zeros((self.M,self.M))

    def solve(self):
        self.timeStart = time()
        while (self.iter < self.maxIters) and (self.err > self.errorTol):
            self._update_field()
            self._write_iteration()
            self.iter += 1
        self.timeEnd = time()
        self.simTime = self.timeEnd - self.timeStart
        self.lambVec.append(self.lamb)
        self.orderParamVec.append(self._calc_orderParam())
        self._writeParams()
        self._write_density()

    def _update_field(self):
        self._fill_F()
        self._calc_error()
        self._update_mixing()
        self._anderson()
        self._get_X()
    
    def _update_mixing(self):
        if self.iter > 200 and self.err < 1e-2:
            self.mixCoeff = 5e-5
        else:
            self.mixCoeff = self.mix0

    def _fill_X(self):
        self.X[0:self.M**2] = self.phi.flatten()
        if self.continuation:
            self._set_continuation()

    def _get_X(self):
        self.phi = self.X[0:self.M**2].reshape((self.M,self.M))
        if self.continuation:
            self._get_continuation()

    def _fill_F(self):
        for i in range(1,self.M-1):
            for j in range(1,self.M-1):
                self.F2d[i,j] = (self.phi[i,j+1] - 2*self.phi[i,j] + self.phi[i,j-1])/self.h**2 + (self.phi[i+1,j] - 2*self.phi[i,j] + self.phi[i-1,j])/self.h**2 + self.lamb*np.exp(self.phi[i,j])

        # Left boundary
        self.F2d[0,:] = (self.phi[1,:] - 2*self.phi[0,:])/self.h**2 + self.lamb*np.exp(self.phi[0,:])

        # Right boundary
        self.F2d[-1,:] = (self.phi[-2,:] - 2*self.phi[-1,:])/self.h**2 + self.lamb*np.exp(self.phi[-1,:])

        # Upper boundary
        self.F2d[:,0] = (self.phi[:,1] - 2*self.phi[:,0])/self.h**2 + self.lamb*np.exp(self.phi[:,0])

        # Lower boundary
        self.F2d[:,-1] = (self.phi[:,-2] - 2*self.phi[:,-1])/self.h**2 + self.lamb*np.exp(self.phi[:,-1])

        # Continuation entry
        if self.continuation:
            self.F[:-1] = self.F2d.flatten()
            if self.stepNum <= self.numKeep:
                self.F[-1] = 0.0
            elif self.continuation:
                self.F[-1] = self.dS - np.dot(self.dYdS,self.X-self.X0)
        else:
            self.F = self.F2d.flatten()

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

        if (self.anderIter > self.anderM and self.iter > 500 and self.stepNum == 0) or (self.anderIter > self.anderM and self.iter > 250 and self.stepNum > 0):
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

        self.anderIter += 1

        if self.anderIter > 5000:
            self.anderIter = 0

    def step(self):
        # Store values from step
        self.iter = 0
        self.err = 1
        self.stepNum += 1

        # Take step in parameter
        if not self.continuation:
            self.lamb += self.dlamb
        else:
            self._step_continuation()

    def _calc_orderParam(self): 
        self.orderParamCurrent = np.max(self.phi)
        return self.orderParamCurrent
    
    def _get_continuation(self):
        self.lamb = self.X[-1]
    
    def _set_continuation(self):
        self.X[-1] = self.lamb

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
            if np.abs(self.dYdS[-1]) > 0.1:
                self.dS *= 1.2
            else:
                self.dS *= 0.5

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
                indices = np.array([prevIndx,prevIndx2,prevIndx3])
                self.lambdac = np.max(self.X_Store[-1,indices])

        # Store continuation parameter
        self._get_continuation()