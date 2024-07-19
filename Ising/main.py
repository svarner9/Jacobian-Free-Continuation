#============================================
#  Interacting Ising Model with External Field
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
from ising import Ising

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, help='Input file for Ising',default='input.json',required=False,dest='inputFile')
args = parser.parse_args()
inputFile = args.inputFile

ising = Ising(inputFile)

startTime = time.time()
for i in tqdm(range(ising.steps)):
    ising.solve()
    ising.step()
print('Total Program Time :', round(time.time() - startTime,3))