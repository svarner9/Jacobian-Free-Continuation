#============================================
#  2D Bratu Problem
#       Solved using Jacobian-free PAC
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

from bratu import Bratu

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help='Input file for Bratu',default='input.json',required=False,dest='inputFile')
args = parser.parse_args()
inputFile = args.inputFile

solver = Bratu(inputFile)

startTime = time.time()
print('Starting calculation...')
for i in tqdm(range(solver.steps)):
    solver.solve()
    solver.step()
print('Total Program Time :', round(time.time() - startTime,3))
