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
from pathlib import Path

# Angstrom = 1/0.5291772
# bins = 200
# z_edges = np.linspace(-15., 15., bins+1)
# z_mid = 0.5*(z_edges[:-1] + z_edges[1:])
# dz = z_edges[1] - z_edges[0]

for potential_type in ['NNP','NNPext']:

    os.chdir(potential_type)

    for H5Fname in ['test-U-5.0.h','test-U+5.0.h']:


        # fs = sorted(glob('*jdftxout'))  # remember to change here !!!
        # fs = sorted(glob('md.o10782753'))  
        # fs = sorted(glob('md.o1078xxxxx'))  
        # for path in Path('src').rglob('*.c'):
            # print(path.name)

        # print('rglob',Path().rglob('*-5.0.h'))

        # fs = sorted(glob('*/*-5.0.h'))  # change first equilibration period to not have this name
        # print("glob",fs)
        # fs = sorted(glob('NaClplanarNaU500Sig2.jdftxout'))  # remember to change here !!!
        # os.path.basename(glob('./Run4/*/NaCl*')[0])
        AllDensities = []  # collect densities from each file after equilibration

        for iFile, f in enumerate(Path().rglob(H5Fname)):
            # print(f)
            # sys.exit(0)
            with h5py.File(f, "r") as fp:
                r = np.array(fp["r"])
                n = np.array(fp["n"])
            AllDensities.append(n.copy())

        AllDensities = np.array(AllDensities)
        AllDensitiesOrig = AllDensities
        # print(np.shape(AllDensities))
        # sys.exit(0)

        samples = iFile + 1 #last index of file (0 base) + 1
        AllDensitiesMean = np.mean(AllDensitiesOrig, axis=0)
        AllDensitiesStd = np.std(AllDensitiesOrig, axis=0)/np.sqrt(samples) # standard error
        # print(np.shape(AllDensitiesMean))

        # save to HDF file


        with h5py.File('All'+H5Fname, "w") as fp:
                        fp["r"] = r
                        fp["n"] = AllDensitiesMean
                        fp["std"] = AllDensitiesStd
                        # fp["V"] = self.force_callback.get_potential()
                        # fp.attrs["T"] = 1300
                        # if self.P is not None:
                        # fp.attrs["P"] = -1 # none doesn't work
                        # fp.attrs["geometry"] = self.force_callback.geometry_type.__name__

    os.chdir('..')

# baseName = os.path.basename(inFile)[-2:]

# fig6 = plt.figure(6)
# fig6.clear(True)
# for i,atName in enumerate(np.unique(atNames)):
#     # average densities over all timesteps
#     # could Reject first 500 steps for equilibration :500
#     # density = np.mean(densityOrig[:,i], axis=0)  # individual file density
#     density = AllDensitiesMean[:,i]  # all file density
#     err = AllDensitiesStd[:,i]
#     # print(len(density))
#     # print(z_mid[-int(bins/2):])
#     # density = gaussian_filter1d(density, 3)
#     x = z_mid[-int(bins/2):]
#     plt.plot(x, density, label=atName)
#     plt.fill_between(x, density-err, density+err,facecolor='r',alpha=0.5)

#     plt.xlabel("z")
#     plt.ylabel("Density")        
#     # plt.xlim((-0.5,0.5))
#     # plt.ylim((0.,4))
#     plt.legend()
#     # fig6.savefig('Plot_NaCl_dens_out_{}.png'.format(baseName))
#     fig6.savefig('Plot_NaCl_dens_out_ALL.png')

