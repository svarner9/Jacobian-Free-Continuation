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

import sys
import time
import argparse
from tqdm import tqdm

sys.path.append('src')
from gsd import GroundStateDominance

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help='Input file for Ising',default='input.json',required=False,dest='inputFile')
args = parser.parse_args()
inputFile = args.inputFile

model = GroundStateDominance(inputFile)

startTime = time.time()
print('Starting calculation...')
for i in tqdm(range(model.steps)):
    model.solve()
    model.step()
print('Total Program Time :', round(time.time() - startTime,3))