#============================================
#  2D Bratu Problem
#       Solved using PAC with Newton
#
#   Written by Sam Varner and Chris Balzer (2024)
#   
#   See top-level README.md for citation.
#   Distributed under the MIT License.
#============================================

import sys
import argparse

sys.path.append('src')
from newton import NewtonBratu

parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Gridpoints in 1-dimension (NxN mesh)',default=50,required=False,dest='N')
parser.add_argument('-dS', type=float, help='Step size in pseudo-arclength continuation',default=0.5,required=False,dest='dS')
args = parser.parse_args()
N = args.N
dS = args.dS

# Initialize class and run
engine = NewtonBratu(N = N, dS = dS)
engine.run()