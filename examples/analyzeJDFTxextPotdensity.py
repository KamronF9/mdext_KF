#!/usr/bin/env python
# read in jdftxout file and produce a density profile plot 
# module load python
# python /global/homes/k/kamron/mdext_KF/examples/analyzeJDFTxextPotdensity.py

import numpy as np
import sys
import gzip
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from glob import glob
import os

Angstrom = 1/0.5291772
z_edges = np.linspace(-10., 10., 201)
z_mid = 0.5*(z_edges[:-1] + z_edges[1:])
dz = z_edges[1] - z_edges[0]

# fs = sorted(glob('*jdftxout'))  # remember to change here !!!
fs = sorted(glob('md.o10782753'))  # remember to change here !!!
# fs = sorted(glob('md.o1078xxxxx'))  # remember to change here !!!
# fs = sorted(glob('md.o*'))  # remember to change here !!!
# fs = sorted(glob('NaClplanarNaU500Sig2.jdftxout'))  # remember to change here !!!
# os.path.basename(glob('./Run4/*/NaCl*')[0])
for inFile in fs:
# for inFile in sys.argv[1:]:
    #Read JDFTx input and convert to VASP
    nSteps = 0 #number of processed steps
    latvecActive = False #Whether reading lattice vectors
    stressActive = False #Whether reading stress tensor
    atposActive = False #Whether reading atomic positions
    PE = []
    T = []
    P = []
    vol = []
    density = []
    fpIn = gzip.open(inFile,'rt') if (inFile[-3:]==".gz") else open(inFile,'r')
    for iLine,line in enumerate(fpIn):
        if line.find('total atoms') > 0:
            nAtoms = int(line.split()[4])
        if line.startswith('IonicDynamics: Step:'):
            tokens = line.split()
            PE.append(float(tokens[4]))
            T.append(float(tokens[8]))
            P.append(float(tokens[10]))
        #Lattice vectors:
        if latvecActive and iLine<refLine+3:
            iRow = iLine-refLine
            R[iRow] = [ float(tok)/Angstrom for tok in line.split()[1:-1] ]
            if iRow==2:
                latvecActive = False
                vol.append(np.abs(np.linalg.det(R)))
        if line.startswith('R ='):
            latvecActive = True
            refLine = iLine+1
            R = np.zeros((3,3))

        #Atomic positions:
        if atposActive and iLine<refLine+nAtoms:
            iRow = iLine-refLine
            tokens = line.split()
            atNames.append(tokens[1])
            atpos[iRow] = [ float(tok) for tok in tokens[2:5] ]
            if iRow+1==nAtoms:
                atposActive = False
                # assert coordsType != "cartesian"
                atposOrig = atpos
                # reset atpos=atposOrig
                if coordsType == "cartesian":
                    # scaled cartesian to fractional coords
                    # convert to fractional using Bohrs (not Ang) OR put atpos move to ang
                    atpos = np.dot(atpos/Angstrom,np.linalg.inv(R.T)) 
                atpos -= np.floor(0.5 + atpos) #wrap
                atpos = np.dot(atpos, R.T) #convert to Cartesian (Angstrom)
                # np.dot(atpos[0],np.linalg.inv(R.T))
                # np.matmul(atpos[0], R.T) # same
                atNames = np.array(atNames)
                # calculate density for each atname - assume in order?
                # for atName in np.unique(atNames):
                # for i in np.argsort(atNames):
                #   print(i)
                # for atName in atNames:
                atUniqueInd = []
                AtUniDensity = []
                for atName in np.unique(atNames):
                    atUniqueInd.append(np.where(atNames == atName))
                    
                for i in range(len(np.unique(atNames))):
                    # print(i)
                    AtUniDensity.append(np.histogram(atpos[atUniqueInd[i], 2], z_edges)[0] / dz)
                    # density.append(np.histogram(atpos[:, 2], z_edges)[0] / dz)  
                    # focus on the z dimension where the potential is applied along xy evenly
                density.append(np.stack(AtUniDensity))
                nSteps += 1
                
        if line.startswith('# Ionic positions in '):
            atposActive = True
            refLine = iLine+1
            atpos = np.zeros((nAtoms,3))
            atNames = []
            coordsType = line.split()[4]
            
    density = np.array(density)
    densityOrig = density
    density = densityOrig
    
    baseName = os.path.basename(inFile)[-2:]
  
    fig6 = plt.figure(6)
    fig6.clear(True)
    for i,atName in enumerate(np.unique(atNames)):

        density = np.mean(densityOrig[:,i], axis=0)  #Reject first 500 steps for equilibration :500
        density = gaussian_filter1d(density, 3) 
        plt.plot(z_mid, density, label=atName)
        plt.xlabel("z")
        plt.ylabel("Density")        
        # plt.xlim((-0.5,0.5))
        # plt.ylim((0.,4))
        plt.legend()
        fig6.savefig('Plot_NaCl_dens_out_{}.png'.format(baseName))
        

