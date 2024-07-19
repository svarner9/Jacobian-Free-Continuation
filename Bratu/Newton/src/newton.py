#============================================
#  2D Bratu Problem
#       Solved using PAC with Newton
#
#   Written by Sam Varner and Chris Balzer (2024)
#   
#   See top-level README.md for citation.
#   Distributed under the MIT License.
#============================================

import os
import time
import numpy as np
from scipy.sparse import lil_matrix

class NewtonBratu:
    def __init__(self,N,dS=1.0):
        # Set parameters
        self.N           = N
        self.h           = 1/N
        self.s           = 0
        self.dS          = dS
        self.steps       = 1000

        # Set initial point
        self.x0          = np.zeros((N**2 + 1,1))

        # Convenient to precompute
        self.lapl        = self.create_laplacian_matrix()
        self.eye         = np.eye(self.N**2, self.N**2)

        # Output
        self.outputDir = 'output/'
        self.output      = np.zeros((self.steps+1,3))
        self.output[0,:] = [self.s, self.x0[-1][0], self.measure(self.x0)]
        self.turnCount = 0
        if not os.path.exists(self.outputDir):
            os.makedirs(self.outputDir)

        print("==================================================")
        print("Bratu Problem in 2D")
        print("  Δu + λ exp(u) = 0")
        print("  Pseudo-arclength continuation with Newton solver")
        print("==================================================")
        print("  Grid size    = (%.f x %.f)"%(self.N,self.N))
        print("  Matrix sizes = (%.f x %.f)"%(self.N**2,self.N**2))
        print("  Contour step = %.2e"%(self.dS))

    def run(self):
        print("==================================================")
        print("Running.....")
        print("==================================================")
        
        avg_iters = 0
        start_time = time.time()
    
        # Take steps
        for self.stepNum in range(1,self.steps+1):
            # Calculate tangent and solve step
            if self.stepNum == 1:
                self.xDot = self.tangent(self.x0)
            else:
                self.xDot = self.tangent(self.x0,old=self.xDot)
            self.x, newton_iters  = self.predictor_corrector(self.x0,self.xDot,self.dS)
            
            # Prepare for next iteration
            self.x0 = np.copy(self.x)

            # Newton iteration counter
            avg_iters +=  newton_iters

            # Store
            self.s += self.dS
            self.output[self.stepNum,:] = [self.s, self.x[-1][0], self.measure(self.x)]

            if self.stepNum > 2:
                if np.sign(self.output[self.stepNum,1]-self.output[self.stepNum-1,1]) != np.sign(self.output[self.stepNum-1,1]-self.output[self.stepNum-2,1]):
                    self.turnCount += 1

            # Printing
            if self.stepNum % 10 == 0:
                print("    Step #: %.f, lambda: %.4f, norm: %.4f"%(self.stepNum,self.output[self.stepNum,1],self.output[self.stepNum,2]))
            
            self.writeParams()

            # Check max temperature
            if self.output[self.stepNum,2] > 10:
                self.output = self.output[:self.stepNum+1]
                break
        
        print("Done!")
        print("  Total time        ---> %.4f seconds"%(time.time() - start_time))
        print("  Total steps       ---> %.f"%(self.stepNum))
        print("  Avg. Newton Steps ---> %.2f per iteration"%(avg_iters/self.stepNum))
        print("  Critical lambda   ---> %.4f"%np.max(self.output[:,1]))
        print("==================================================")

    def measure(self,x):
        return np.max(x[:-1])
    
    def create_laplacian_matrix(self):
        laplacian = lil_matrix((self.N**2,self.N**2))
        for i in range(self.N):
            for j in range(self.N):
                index = i*self.N + j
                if i > 0:
                    laplacian[index, index - self.N] = 1
                if i < self.N-1:
                    laplacian[index, index + self.N] = 1
                if j > 0:
                    laplacian[index, index-1] = 1
                if j < self.N-1:  
                    laplacian[index, index+1] = 1
                laplacian[index, index] = -4
        laplacian /= self.h**2
        
        # Apply boundary condition
        for i in range(self.N):
            for j in range(self.N):
                index = i*self.N + j
                if i == 0 or i == self.N-1 or j == 0 or j == self.N-1:
                    laplacian[index,:] = 0
                    laplacian[:,index] = 0
                    laplacian[index,index] = 1
        return laplacian.tocsr()
    
    def tangent(self,x,old=None):
        Ju, Jlam = self.jacobian(x)
        xDot = -np.linalg.solve(Ju,Jlam)
        xDot = np.append(xDot,1)
        xDot /= np.linalg.norm(xDot)
        xDot = xDot[:,None]
        if not old is None:
            if np.dot(np.transpose(xDot),old) < 0:
                xDot *= -1
        return xDot
    
    def residual(self,x,x0,xDot,ds):
        Fu    = self.lapl*x[:-1] + x[-1]*np.exp(x[:-1])  # size N^2 x 1
        Flam  = np.dot(np.transpose(xDot),x-x0) - ds     # size 1x1 

         # Adjust at boundary
        for i in range(self.N):
            for j in range(self.N):
                index = i*self.N + j
                if i == 0 or i == self.N-1 or j == 0 or j == self.N-1:
                    Fu[index] = x[index]
        return Fu,Flam
    
    def jacobian(self,x):
        Ju   = self.lapl + x[-1]*np.exp(x[:-1])*self.eye   # size N^2 x N^2
        Jlam = np.exp(x[:-1])                              # size N^2 x 1

        # Adjust at boundary
        for i in range(self.N):
            for j in range(self.N):
                index = i*self.N + j
                if i == 0 or i == self.N-1 or j == 0 or j == self.N-1:
                    Ju[index,index] = 1
                    Jlam[index]     = 0
        return Ju,Jlam
    
    def predictor_corrector(self,x0,xDot0,ds,tol=1e-6,max_iter=100):
        x = x0 + xDot0*ds
        for iteration in range(max_iter):
            Ju, Jlam = self.jacobian(x)
            Fu, Flam = self.residual(x,x0,xDot0,ds)
            Jinv = np.linalg.inv(np.concatenate((np.concatenate((Ju, Jlam),axis=1), np.transpose(xDot0))))
            step =  -Jinv*np.concatenate((Fu,Flam),axis=0)
            x   += step

            if np.linalg.norm(step) < tol:
                break
        if iteration == max_iter-1:
            print("Newton did not converge")
        return x, iteration

    def writeParams(self):
        if self.stepNum == 1:
            tag = 'w'
        else:
            tag = 'a'
        with open(self.outputDir + 'params.dat',tag) as f:
            f.write(str(self.stepNum) + ' ' + f'{self.output[self.stepNum,1]:.5e}' + ' ' + f'{self.output[self.stepNum,2]:.5e}')
            f.write(' ' + f'{self.s:.5f}' + ' ' + f'{self.turnCount:.0f}')
            f.write(' \n')