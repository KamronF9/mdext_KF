#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from io import StringIO
import sys
import os

# print('units in ang need change')
print('assume only xy tilt and 0 lo bounds for now')
Angstrom = 1/0.5291772
# Angstrom = 1  # HACK for validation

if len(sys.argv) < 6:
    print('Usage: lammpsToJDFTx.py <dumpFile> <outPrefix> <step stride> <total snapshots take> <current run folder index number>')
    # lammpsToJDFTx.py round11extraHr3_8kOOljOn.dump test 1 2 
    # python /home/kamron/mdext_KF/examples/lammpsToJDFTxTricOrOrtho.py pylammps+5.0last.dump init+5.0 1 1 1 
    # python /home/kamron/mdext_KF/examples/lammpsToJDFTxTricOrOrtho.py pylammps-5.0last.dump init-5.0 1 1 1 
    # python /home/kamron/mdext_KF/examples/lammpsToJDFTxTricOrOrtho.py pylammps+2.5last.dump init+2.5 1 1 1 
    # python /home/kamron/mdext_KF/examples/lammpsToJDFTxTricOrOrtho.py pylammps-2.5last.dump init-2.5 1 1 1 
    #/global/homes/k/kamron/mdext_KF/examples/lammpsToJDFTxTricOrOrtho.py
    exit(1)

dumpFile = sys.argv[1]  #read lammps dumpFile
outPrefix = sys.argv[2] #create outPrefix.lattice and outPrefix.ionpos
nEvery = int(sys.argv[3]) #length of stride to take for gathering snapshots to run 
snapsToTake = int(sys.argv[4]) #total snapshots to take
currFolderIndex = int(sys.argv[5]) #current index number for appending to folder - not used here
# now get atom name from dump file
#atomNames = ['Mg','Cl'] #atom names in order of lammps ID # TS uses alphabetical order for mixture**

#Read dump file:
tStepSection = False
nAtomsSection = False
boxSection = False
atomsSection = False
headerWritten = False
stepActive = False #Whether to process current data
writeOut = False # when to writeout to file
triclinic = False #if dump is triclinic then True - ortho if false
snapsTaken = 0

# make run# dir for each ionpos folder
# for i in range(snapsToTake):
#     os.system(f'mkdir run{(currFolderIndex + i):04d}')
#     # print(f'mkdir run{(0+1):04d}')
for line in open(dumpFile):
    #Atom data:
    if atomsSection:
        # atomData[atomLine] = np.loadtxt(StringIO(line))  # id type element x y z fx fy fz
        atomData.append(line.split())  # id type element x y z fx fy fz
        atomLine += 1
        if atomLine == nAtoms:
            atomsSection = False #Frame complete
            writeOut = True
    if line.startswith("ITEM: ATOMS"):
        atomsSection = True
        atomLine = 0
        # atomData = np.zeros((nAtoms,#))
        atomData = []
    #Box:
    if boxSection:
        tokens = line.split()
        # ITEM: BOX BOUNDS xy xz yz pp pp pp
        # 0.0000000000000000e+00 1.2446608774798516e+01 4.1488625359178286e+00
        # 0.0000000000000000e+00 7.1860530340313176e+00 0.0000000000000000e+00
        # 0.0000000000000000e+00 3.0000000000000000e+01 0.0000000000000000e+00
        # xlo_bound = xlo + MIN(0.0,xy,xz,xy+xz)
        # xhi_bound = xhi + MAX(0.0,xy,xz,xy+xz)
        # ylo_bound = ylo + MIN(0.0,yz)
        # yhi_bound = yhi + MAX(0.0,yz)
        # zlo_bound = zlo
        # zhi_bound = zhi
        # HACK xxxxxxx add back in *Angstrom
        
        
        
        if triclinic:
            # Triclinic
            tricBox[boxLine] = [ float(tok) for tok in line.split() ]  # convert to bohr from ang later in  output
            boxLine += 1

            if boxLine == 3:
                boxSection = False
                R = np.zeros((3,3))
                # assume only xy tilt and 0 lo bounds for now
                R[0,0] = tricBox[0,1] - tricBox[0,2]
                R[1,0] = tricBox[0,2]
                R[1,1] = tricBox[1,1] 
                R[2,2] = tricBox[2,1]
                # print(R)
            
        else:
            # Orthorhombic 
            L[boxLine] = float(tokens[1]) - float(tokens[0])
            boxLine += 1
            if boxLine == 3:
                boxSection = False
                R = np.zeros((3,3))
                # assume only xy tilt and 0 lo bounds for now
                R[0,0] = L[0]
                R[1,1] = L[1]
                R[2,2] = L[2]
                # print(R)
    
    if line.startswith("ITEM: BOX BOUNDS"):
        boxSection = True
        boxLine = 0
        if line.find('xy xz yz') > 0:
            #triclinic
            tricBox = np.zeros((3,3))
        else:
            # orthorhombic
            L = np.zeros(3)
    #Number of atoms:
    if nAtomsSection:
        nAtoms = int(line)
        nAtomsSection = False
    if line.startswith("ITEM: NUMBER OF ATOMS"):
        nAtomsSection = True
    #Time step:
    if tStepSection:
        iStep = int(line)
        # stepActive = (iStep % nEvery == 0 and not iStep == 0)
        stepActive = (iStep % nEvery == 0)
        tStepSection = False
    if line.startswith("ITEM: TIMESTEP"):
        tStepSection = True


    
    #Write lattice and atoms in JDFTx format:
    if stepActive and writeOut and snapsTaken<snapsToTake:
        print(iStep)
        stepActive = False
        writeOut = False
        
        atomData = np.array(atomData)
    # print(np.array(atomData).T)
    # print(np.lexsort((np.array(atomData).T)[:2]))

        #transpose and convert to Bohrs
        R = R.T * Angstrom

        with open(outPrefix+'.lattice', 'w') as fp:
            print('lattice \\', file=fp)
            print('{:15.8f}{:15.8f}{:15.8f} \\'.format(*tuple(R[0])), file=fp)
            print('{:15.8f}{:15.8f}{:15.8f} \\'.format(*tuple(R[1])), file=fp)
            print('{:15.8f}{:15.8f}{:15.8f}'.format(*tuple(R[2])), file=fp)
                # lattice \
                # 	 23.52068 0 0 \
                # 	0  13.57967 0 \
                # 	0 0  56.69178
            # fp.write(f'lattice \\\n\t{:9.5f} 0 0 \\\n\t0 {:9.5f} 0 \\\n\t0 0 {:9.5f}\n'.format(*tuple(R)))
        atomData = atomData[np.lexsort((atomData.T)[:2])] #sort by atom type then ID

        with open(outPrefix+'.ionpos', 'w') as fp: # +str(iStep)
            fp.write('# Ionic positions in cartesian coordinates:\n')
            # for iAtom, atType in enumerate(atomData[:,1].astype(int)-1):
            for iAtom in range(nAtoms):
                fp.write('ion {:s} {:9.5f} {:9.5f} {:9.5f} 1\n'.format(atomData[iAtom,2], *tuple(atomData[iAtom,3:6].astype(float)*Angstrom)))
                # fp.write('ion {:s} {:s} {:9.5f} {:9.5f} {:9.5f} 1\n'.format(atomData[iAtom,0], atomData[iAtom,2], *tuple(atomData[iAtom,3:6].astype(float)*Angstrom)))
        
        # move into each run folder
        # os.system(f'mv -t run{currFolderIndex + snapsTaken}/ {outPrefix}.lattice {outPrefix}.ionpos')
        snapsTaken += 1