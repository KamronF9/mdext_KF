#!/usr/bin/env python
# To read in jdftxout file and produce a density profile plot 
#
# module load python
# python /global/homes/k/kamron/mdext_KF/examples/analyzeJDFTxextPotdensity.py

import numpy as np
import sys
import gzip
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from glob import glob
import os
import h5py

Angstrom = 1/0.5291772
bins = 200
z_edges = np.linspace(-15., 15., bins+1)
z_mid = 0.5*(z_edges[:-1] + z_edges[1:])
dz = z_edges[1] - z_edges[0]


H5Fname = 'test-U+5.0.h'
# H5Fname = 'test-U-5.0.h'

# fs = sorted(glob('*jdftxout'))  # remember to change here !!!
# fs = sorted(glob('md.o10782753'))  
# fs = sorted(glob('md.o1078xxxxx'))  
fs = sorted(glob('md.o*'))  # change first equilibration period to not have this name
# fs = sorted(glob('NaClplanarNaU500Sig2.jdftxout'))  # remember to change here !!!
# os.path.basename(glob('./Run4/*/NaCl*')[0])
AllDensities = []  # collect densities from each file after equilibration

for iFile,inFile in enumerate(fs):
# for inFile in sys.argv[1:]:
    #Read JDFTx input
    print(f'{inFile=}')
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
                    AtUniDensityTwoSide = np.histogram(atpos[atUniqueInd[i], 2], z_edges)[0] / (dz*R[0,0]*R[1,1])
                    # print(f'{len(AtUniDensityTwoSide)=}') #200
                    # plt.plot(z_mid, AtUniDensityTwoSide)
                    # plt.savefig('test1.png')
                    # plt.plot(z_mid[-int(bins/2):], np.flip(AtUniDensityTwoSide[:int(bins/2)]))
                    # plt.plot(z_mid[-int(bins/2):], AtUniDensityTwoSide[-int(bins/2):])
                    # plt.savefig('test2.png')
                    # sys.exit(0)
                    AtUniDensityOneSide = (np.flip(AtUniDensityTwoSide[:int(bins/2)])+AtUniDensityTwoSide[-int(bins/2):])*0.5
                    # print(f'{len(AtUniDensityOneSide)=}') # 100
                    # print(z_mid[-int(bins/2):])
                    # print(np.flip(z_mid[:int(bins/2)]))
                    # sys.exit(0)
                    AtUniDensity.append(AtUniDensityOneSide)
                    # need to scale the density to the volume to be inline with mdext
                    # AtUniDensity.append(np.histogram(atpos[atUniqueInd[i], 2], z_edges)[0] / dz)
                    # density.append(np.histogram(atpos[:, 2], z_edges)[0] / dz)  
                    # focus on the z dimension (2) where the potential is applied along xy evenly
                # Reflect densities about 0 and average them
                
                # density.append(np.stack(AtUniDensity)) # histogram bin values dimension 2 deep for 2 unique atom types in this case
                AllDensities.append(np.stack(AtUniDensity))
                nSteps += 1
                
        if line.startswith('# Ionic positions in '):
            atposActive = True
            refLine = iLine+1
            atpos = np.zeros((nAtoms,3))
            atNames = []
            coordsType = line.split()[4]
    # if iFile==1: break


# shifting below to be after all processing of files
# density = np.array(density)
# densityOrig = density
# density = densityOrig ?

AllDensities = np.array(AllDensities)
AllDensitiesOrig = AllDensities
# print(np.shape(AllDensities))

samples = iFile + 1 #last index of file (0 base) + 1
AllDensitiesMean = np.mean(AllDensitiesOrig, axis=0).T
AllDensitiesStd = np.std(AllDensitiesOrig, axis=0).T/np.sqrt(samples) # standard error
# print(np.shape(AllDensitiesMean))

# save to HDF file


with h5py.File(H5Fname, "w") as fp:
                fp["r"] = z_mid[-int(bins/2):]
                fp["n"] = AllDensitiesMean
                fp["std"] = AllDensitiesStd
                # fp["V"] = self.force_callback.get_potential()
                fp.attrs["T"] = 1300
                # if self.P is not None:
                fp.attrs["P"] = -1 # none doesn't work
                # fp.attrs["geometry"] = self.force_callback.geometry_type.__name__

# baseName = os.path.basename(inFile)[-2:]

fig6 = plt.figure(6)
fig6.clear(True)
for i,atName in enumerate(np.unique(atNames)):
    # average densities over all timesteps
    # could Reject first 500 steps for equilibration :500
    # density = np.mean(densityOrig[:,i], axis=0)  # individual file density
    density = AllDensitiesMean[:,i]  # all file density
    err = AllDensitiesStd[:,i]
    # print(len(density))
    # print(z_mid[-int(bins/2):])
    # density = gaussian_filter1d(density, 3)
    x = z_mid[-int(bins/2):]
    plt.plot(x, density, label=atName)
    plt.fill_between(x, density-err, density+err,facecolor='r',alpha=0.5)

    plt.xlabel("z")
    plt.ylabel("Density")        
    # plt.xlim((-0.5,0.5))
    # plt.ylim((0.,4))
    plt.legend()
    # fig6.savefig('Plot_NaCl_dens_out_{}.png'.format(baseName))
    fig6.savefig('Plot_NaCl_dens_out_ALL.png')

