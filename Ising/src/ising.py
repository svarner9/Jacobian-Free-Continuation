#============================================
#  Interacting Ising Model with External Field
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

class Ising(object):
    def __init__(self, inputFile):

        # Default input file parameters
        self.a = 3.0
        self.h0 = -2.0
        self.dh = 0.01
        self.steps = 5000
        self.anderM  = 3
        self.mixCoeff = 1e-5
        self.errorTol = 1.0e-6
        self.continuation = True
        self.outputDir = './'
        
        # Read input file
        self.inputFile = 'input/' + inputFile
        self._read_input()

        # Parameters
        self.m_vec = np.zeros((self.steps))
        self.m = 0.9999999 * np.sign(self.h0)

        # Stepping
        self.stepNum = 0
        self.h = self.h0
        self.h_vec = np.zeros((self.steps))
        
        # Mixing
        self.iter = 0
        self.err = 1.0
        self.maxIters = 1e7

        # Arrays for solvers
        self.numKeep = 10
        if self.continuation:
            self.numUnknowns = 2
        else:
            self.numUnknowns = 1
        self.X        = np.zeros((self.numUnknowns,))
        self.F        = np.zeros((self.numUnknowns,))        
        self.X_Store  = np.zeros((self.numUnknowns,self.numKeep))
        self.F_Store  = np.zeros((self.numUnknowns,self.numKeep))

        # Anderson variables
        self.X_And  = np.zeros((self.numUnknowns,self.anderM))
        self.G_And  = np.zeros((self.numUnknowns,self.anderM))

        # Continuation variables
        self.s           = 0.0
        self.dS          = self.dh
        self.order       = 0.0
        self.turnCount   = 0
        self.s_vec       = np.zeros((self.steps))
        self.orderParam  = np.zeros((self.steps))
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

    def writeParams(self):
        if self.stepNum == 0:
            tag = 'w'
        else:
            tag = 'a'
        with open(self.outputDir + 'params_.dat',tag) as f:
            f.write(str(self.stepNum) + ' ' + f'{self.h:.5f}' + ' ' + f'{self.m:.5f}')
            if self.continuation:
                f.write(' ' + f'{self.s:.5f}' + ' ' + f'{self.turnCount:.0f}')
            f.write(' \n')

    def solve(self):
        while (self.iter < self.maxIters) and (self.err > self.errorTol):
            self._update_field()
            self.iter += 1
        self.writeParams()

    def _update_field(self):
        self._fill_F()
        self._calc_error()
        self._anderson()
        self._get_X()

    def _fill_X(self):
        self.X[0] = self.m
        if self.continuation:
            self._set_continuation()

    def _get_X(self):
        self.m = self.X[0]
        if self.continuation:
            self._get_continuation()

    def _fill_F(self):
        self.F[0] = -(-self.a * self.m - self.h + 0.5 * np.log((1+self.m)/(1-self.m)))
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
                print(' Lstsq Failed')
                self._fill_X()                
                self._picard()
        else:
            self._picard()

    def step(self):
        # Store values from step
        self.h_vec[self.stepNum] = self.h
        self.m_vec[self.stepNum] = self.m
        
        self.iter = 0
        self.err = 1
        self.stepNum += 1

        # Take step in parameter
        if not self.continuation:
            self.h += self.dh
        else:
            self._step_continuation()

    def _get_continuation(self):
         self.h = self.X[-1]
    
    def _set_continuation(self):
        self.X[-1] = self.h

    def _calc_orderParam(self):
        return self.m

    def _step_continuation(self):
        # Convenient indices
        currentIndx = self.stepNum % self.numKeep
        prevIndx    = (self.stepNum - 1) % self.numKeep
        prevIndx2   = (self.stepNum - 2) % self.numKeep
        prevIndx3   = (self.stepNum - 3) % self.numKeep

        # Store previous profile
        self.s_vec[prevIndx]       = self.s
        self.orderParam[prevIndx]  = self._calc_orderParam()
        self.X_Store[:,prevIndx]   = self.X
        self.X0                    = deepcopy(self.X)

        # Take predictor step
        if self.stepNum <= self.numKeep:
            # Natural continuation step
            self.X[-1]  += self.dS
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
            minDs = 0.01
            maxDs = 0.1
            if self.dS < minDs:
                self.dS = minDs
            elif self.dS > maxDs:
                self.dS = maxDs

            # Take step
            self.X  += self.dYdS*self.dS

            # Increment s
            self.s   += self.dS

           # Determine turn count
            if (self.X_Store[-1,prevIndx] - self.X_Store[-1,prevIndx2])/ds1  * (self.X_Store[-1,prevIndx2] - self.X_Store[-1,prevIndx3])/ds2 < 0:
                self.turnCount += 1

        # Store continuation parameter
        self._get_continuation()