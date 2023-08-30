import numpy as np
import pandas as pd
import argparse
import glob
import os
import sys

curr_dir = os.getcwd()

parser = argparse.ArgumentParser(description="write lammps thermo output to a csv file")
parser.add_argument('-i', action="store", dest="input_file")
args = parser.parse_args()

# purge all empty strings from token lists - and return list
def purge(tokens):
        return [t for t in tokens if len(t) >= 1]

def is_num(s):
    try:
        float(s)
    except ValueError:
        return False
    else:
        return True

# Step Temp PotEng TotEng Press Volume lambda lj dlj 1 gT qNa qCl
# Step Temp PotEng TotEng Press Volume lambda lj dlj 1

# Step Temp PotEng TotEng Press Volume lambda lj dlj 1 qNa qCl

# ======== READ INPUT FILE ===========

infile = args.input_file

with open(infile,'r') as f:
        lines = f.readlines()

lines = [l.strip('\n') for l in lines]

datastart_idxs = []
dataend_idxs = []

linecounter = 0
for line in lines:
        # if line.startswith('Step'):
        if line.startswith('   Step'):
                datastart_idxs.append(linecounter)
        linecounter += 1
# print(f'{datastart_idxs=}')
linecounter = 0
for line in lines:
        if (line.startswith('Loop time of')) and (lines.index(line) > np.min(datastart_idxs)):
                dataend_idxs.append(linecounter)
        linecounter += 1

print(datastart_idxs)
print(dataend_idxs)

# may not need with the reset with mdext.... 
# print(purge(lines[datastart_idxs[0]].split(' ')))
# sys.exit(0)
header = purge(lines[datastart_idxs[1]].split(' ')) # should be the same header for all  but isn't!! so use second one
# we only want the TI data not initial equilibration
data = []
for idx in range(len(datastart_idxs)): # index of indices
        data_chunk = lines[datastart_idxs[idx]+1:dataend_idxs[idx]]
        for d in data_chunk:
                # print(purge(d.split(' ')))
                # sys.exit(0)
                data.append(purge(d.split(' ')))

df = pd.DataFrame(data, columns=header)
outfile = infile[:-6] + '.csv'
df.to_csv(outfile)