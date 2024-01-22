#!/usr/bin/env python
# RDF only
import numpy as np
import sys
import h5py
import glob
import matplotlib.pyplot as plt
# import mldft1d

# read in H5 files for each ext potential

# 1D LJ:
# endRange = 15.0
# stepSize = 0.5
# decimals = 1

# Water:
# Define particle if multiple particles used - H2O
particle = 1 # 0 based so 1 is O in H2O
endRange = 0.1
stepSize = 0.01
decimals = 2

filename = 'AllData003water.h'

rAll = []
nAll = []
VAll = []
EAll = []
lbdaAll = []

def mirrorData(data,sign):
    lenData = len(data)
    mirroredData=np.zeros(2*lenData-1)
    # print(sign*data[::-1][:-1])
    mirroredData[:lenData-1]=sign*data[::-1][:-1] # reverse and leave out last (ie 0 pt)
    mirroredData[lenData-1:]=data

    # n_sym=np.zeros((2*len(n),2))
    # np.zeros_like(n)
    # n_sym[:len(n)]=n[::-1]
    # n_sym[len(n):]=n
    return mirroredData

for Ui in np.around(np.arange(0,endRange*2 + stepSize ,stepSize), decimals=decimals)-endRange:  
    # print(f"{Ui:+.1f}")
    fname = glob.glob(f"*{Ui:+.2f}*h5")[0]
    print('loading file ', fname)

    with h5py.File(fname, "r") as fp:
        r = np.array(fp["r"])
        mir_r = mirrorData(r,-1)

        if 'particle' in locals():
            n = np.array(fp["n"])[:,particle].flatten()
            V = np.array(fp["V"])[:,particle].flatten()
        else:
            n = np.array(fp["n"]).flatten()
            V = np.array(fp["V"]).flatten()
        # print(V)
        # print(np.shape(V)[0])
        # sys.exit(1)
        mir_n = mirrorData(n,1)
        mir_V = mirrorData(V,1)
        # print(mir_r)
        # plt.plot(r,n)
        # plt.savefig('n')
        # plt.clf()
        # plt.plot(mir_r,mir_n)
        # plt.savefig('mir_n')
        # plt.clf()
        # sys.exit(1)

        try:
            E = np.array(fp["pe_history"]).mean()
        except:
            print('no PE data in H5 to read')
            E = np.zeros(1)
        if Ui==0.: 
            print(f'N_bulk is {n.mean()=}')
            n_bulk = n.mean()

    # rAll.append(mir_r.copy())
    nAll.append(mir_n.copy())
    # VAll.append(mir_V.copy())
    EAll.append(E.copy())
    lbdaAll.append(Ui)

# print(np.shape(VAll))
# sys.exit(1)
# generate rolled up H5 file for TI and training

f = h5py.File(filename, "w")
f["z"] = mir_r # was get1D(grid1d.z)
f["V"] = mir_V # will store last shape # was get1D(V.data)  one sample shape
f["lbda"] = lbdaAll # was lbda_arr
f["n"] = nAll  # dim = lambdas, n_points # was n
f["E"] = EAll # was E
# if has_known_part:
#     f["E0"] = E0
#     f["V0"] = V0
# f.attrs["mu"] = dft.mu
f.attrs["n_bulk"] = n_bulk
# for dft_arg_name, dft_arg_value in dft_kwargs.items():
    # f.attrs[dft_arg_name] = dft_arg_value
f.close()



'''
---------------
(mdextKF) [py311] kamron@aimp:~/HardRods1D_KF/examples/data/random_potentials$ h5dump random_n0.4_L1_sigma0.1_seed10.h5 
HDF5 "random_n0.4_L1_sigma0.1_seed10.h5" {
GROUP "/" {
   ATTRIBUTE "R" {
      DATATYPE  H5T_IEEE_F64LE
      DATASPACE  SCALAR
      DATA {
      (0): 0.5
      }
   }
   ATTRIBUTE "T" {
      DATATYPE  H5T_IEEE_F64LE
      DATASPACE  SCALAR
      DATA {
      (0): 1
      }
   }
   ATTRIBUTE "n_bulk" {
      DATATYPE  H5T_IEEE_F64LE
      DATASPACE  SCALAR
      DATA {
      (0): 0.4
      }
   }
   DATASET "E" {
      DATATYPE  H5T_IEEE_F64LE
      DATASPACE  SIMPLE { ( 11 ) / ( 11 ) }
      DATA {
      (0): -0.666667, -0.771767, -0.94304, -1.23197, -1.5696, -1.95642,
      (6): -2.4028, -2.87768, -3.38074, -3.90617, -4.44956
      }
   }
   DATASET "V" {
      DATATYPE  H5T_IEEE_F64LE
      DATASPACE  SIMPLE { ( 20 ) / ( 20 ) }
      DATA {
      (0): -0.125009, -0.239163, -0.571516, -0.808561, -0.727838, -0.300044,
      (6): 0.336844, 0.97065, 1.39311, 1.44471, 1.07465, 0.39691, -0.351899,
      (13): -0.964015, -1.38034, -1.63748, -1.71133, -1.4938, -0.98018,
      (19): -0.412721
      }
   }
   DATASET "lbda" {
      DATATYPE  H5T_IEEE_F64LE
      DATASPACE  SIMPLE { ( 11 ) / ( 11 ) }
      DATA {
      (0): 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5
      }
   }
   DATASET "n" {
      DATATYPE  H5T_IEEE_F64LE
      DATASPACE  SIMPLE { ( 11, 20 ) / ( 11, 20 ) }
      DATA {
      (0,0): 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
      (0,13): 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
      (1,0): 0.364169, 0.384866, 0.451027, 0.504985, 0.486428, 0.397264,
      (1,6): 0.293137, 0.216026, 0.176013, 0.171631, 0.20533, 0.284611,
      (1,12): 0.406459, 0.541034, 0.654334, 0.73514, 0.76123, 0.690012,
      (1,18): 0.543981, 0.416928,
      (2,0): 0.370721, 0.398992, 0.498921, 0.592525, 0.55875, 0.416246,
      (2,6): 0.28202, 0.199, 0.159902, 0.155752, 0.188337, 0.272295,
      (2,12): 0.429873, 0.663611, 0.931942, 1.17777, 1.26774, 1.03307,
      (2,18): 0.670271, 0.445987,
      (3,0): 0.222448, 0.262555, 0.422053, 0.591659, 0.528772, 0.287314,
      (3,6): 0.113974, 0.0448029, 0.0239129, 0.0221867, 0.0383541, 0.104131,
      (3,12): 0.306481, 0.712018, 1.20936, 1.64046, 1.8055, 1.39583, 0.71874,
      (3,19): 0.332292,
      (4,0): 0.156256, 0.195044, 0.367914, 0.577842, 0.49713, 0.219439,
      (4,6): 0.0634352, 0.0181261, 0.00783069, 0.00710432, 0.0147646,
      (4,11): 0.0562737, 0.238794, 0.730363, 1.44132, 2.09029, 2.35778,
      (4,17): 1.72183, 0.73577, 0.266984,
      (5,0): 0.156388, 0.195439, 0.371283, 0.589662, 0.504988, 0.220065,
      (5,6): 0.0633767, 0.0181123, 0.00782698, 0.00710113, 0.0147545,
      (5,11): 0.0562191, 0.239666, 0.753, 1.57112, 2.42457, 2.80837, 1.9249,
      (5,18): 0.759139, 0.268262,
      (6,0): 0.068441, 0.0956876, 0.249918, 0.495573, 0.39375, 0.113737,
      (6,6): 0.017492, 0.00257526, 0.000743816, 0.000650236, 0.0019508,
      (6,11): 0.0144883, 0.129204, 0.68515, 1.78076, 2.8088, 3.2885, 2.23281,
      (6,18): 0.683972, 0.153697,
      (7,0): 0.0436012, 0.0644845, 0.199067, 0.444657, 0.339105, 0.0789216,
      (7,6): 0.00879797, 0.000898399, 0.000240039, 0.000213796, 0.000601035,
      (7,11): 0.00691434, 0.0916869, 0.641941, 1.9041, 3.07268, 3.64742,
      (7,17): 2.42516, 0.63697, 0.112811,
      (8,0): 0.0271689, 0.0427934, 0.156359, 0.393082, 0.287765, 0.0539839,
      (8,6): 0.0042959, 0.00030591, 9.18932e-05, 8.24256e-05, 0.000212341,
      (8,11): 0.00329718, 0.0641329, 0.593653, 2.00724, 3.28425, 3.94213,
      (8,17): 2.58661, 0.585287, 0.0811221,
      (9,0): 0.016893, 0.0279721, 0.121181, 0.343276, 0.24116, 0.0365505,
      (9,6): 0.00205564, 0.000117789, 4.66633e-05, 4.2445e-05, 9.22136e-05,
      (9,11): 0.00162255, 0.0444734, 0.543617, 2.09599, 3.45379, 4.18336,
      (9,17): 2.72327, 0.532245, 0.0582953,
      (10,0): 0.0103221, 0.0179845, 0.0919713, 0.297791, 0.200352, 0.0238541,
      (10,6): 0.000863338, 5.07072e-05, 1.93956e-05, 1.77849e-05,
      (10,10): 3.78016e-05, 0.000858342, 0.0303134, 0.493598, 2.17379,
      (10,15): 3.58998, 4.38076, 2.84137, 0.481046, 0.0414366
      }
   }
   DATASET "z" {
      DATATYPE  H5T_IEEE_F64LE
      DATASPACE  SIMPLE { ( 20 ) / ( 20 ) }
      DATA {
      (0): 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
      (12): 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95
      }
   }
}
}
---------------

# need at least two timesteps
if len(sys.argv)<3:
    print('Usage: RDF.py <lammps dump> <basefile name for RDFs> <frames of start to ignore>')
    exit(1)

# inFile = 'lammps.dump'
# RFile = 'lammps_R.out'
# outFile = 'lammps_proc'
    
inFile = sys.argv[1] #lammps dump file
# RFile = sys.argv[2] #lattice dump file
outFile = sys.argv[2] #basename for RDF files
framesToIgnore = int(sys.argv[3]) # e.g., if 10 and dump saved 1000 then this is 10000th frame

legLabel = '1D LJ'  # label for the RDF plot later

refLine = 0
nSteps = 0 #number of processed steps
nFrames = 0 # number of frames read in
# nEvery = 10 #select this many frames
preprocess = True

latvecActive = False #Whether reading lattice vectors
tricLat = False     #if latvec is orthogonal = False, triclinic =True
stepActive = False #Whether to process current data - tailor later
atomsActive = False #Whether to read total atoms
atposActive = False #Whether reading atomic positions
rdfInited = False #Whether RDF params have been initiated 

# TO DO
# delete first TIMESTEP and moved to end of file

# # For testing lines
# iLine = 0
# f = open(inFile)
# #----
# line = f.readline()
# iLine += 1 
# line

# Preprocess file - to make sure last line has TIMESTEP keyword to know when to stop
if preprocess:
    import subprocess
    def tail(f, n, offset=0):
        proc = subprocess.Popen(['tail', '-n', str(n), f], stdout=subprocess.PIPE)
        lines = proc.stdout.readlines()
        return lines #[:, -offset]
    if not 'TIMESTEP' in str(tail(inFile, 1, 0)):  
        with open(inFile, 'a') as f:  # append mode
            f.write('TIMESTEP')


for iLine,line in enumerate(open(inFile)):
    

    
    #Lattice vectors:
    if latvecActive and iLine<refLine+3:
        iRow = iLine-refLine

        # read each line
        bounds[iRow] = [ float(tok) for tok in line.split() ]
        R = np.array([[bounds[0,1]-bounds[0,0], 0. , 0.],
                        [0., bounds[1,1]-bounds[1,0] , 0.],
                        [0., 0., bounds[2,1]-bounds[2,0] ]])
        if iRow==2:
            latvecActive = False
    if line.startswith('ITEM: BOX BOUNDS'):
        nFrames += 1
        latvecActive = True
        if line.find('xy xz yz') > 0:
            tricLat = True
        else:
            tricLat = False
        refLine = iLine+1
        R = np.zeros((3,3))
        Tric = np.zeros((6))
        bounds = np.zeros((3,2))
    # Atomic positions
    if atposActive and iLine<refLine+nAtoms:
        iRow = iLine-refLine
        tokens = line.split()
        atNames.append(int(tokens[1]))  # index for name either 1 or 2
        atpos[iRow] = [ float(tok) for tok in tokens[2:5] ]  # 2,3,4 
        if iRow+1==nAtoms:
            atposActive = False
            atNames = np.array(atNames)
    if line.startswith('ITEM: ATOMS'):
        atposActive = True
        refLine = iLine+1 # start of where to read in atom positions
        atpos = np.zeros((nAtoms,3))
        atNames = []
    # Number of atoms
    if atomsActive:
        nAtoms = int(line.split()[0])
        atomsActive = False
    if line.find('NUMBER OF ATOMS') > 0:
        atomsActive = True
    # Final processing
    
    if (line.find('TIMESTEP') > 0) and (iLine>5) and (nFrames>framesToIgnore):  # once you get to the end/beginning of the next tally up the RDF
        # RDF initialize
        if not rdfInited:
            rMax = 0.5 * np.mean(np.diag(R))     
            dr = 0.01
            rBins = np.arange(0., rMax, dr)
            rBins[0] = 0.01*dr #ignore self
            rMid = 0.5*(rBins[:-1]+rBins[1:])
            # binVol = (4*np.pi/3)*(rBins[1:]**3 - rBins[:-1]**3)  # if 3D and spherical
            binVol = (rBins[1:] - rBins[:-1]) # if bins are rectangular with x and y = 1
            rdf = np.zeros((len(rMid),3))
            rdfInited = True

        x = np.dot(atpos, np.linalg.inv(R.T))   # normalize positions to lattice shape
        x1 = x[np.where(atNames==1)[0]]  # was Mg/Na
        # print(len(x1))  # number of atoms
        # x2 = x[np.where(atNames==2)[0]]  # was Cl
        def getRDF(x1, x2):
            dx = x1[None,:,:] - x2[:,None,:]  # None adds a dimension 
            dx -= np.floor(0.5+dx) #minimum image convention
            r = np.linalg.norm(np.dot(dx, R.T), axis=-1).flatten()  
            # maybe done to cast relative coords back onto coord basis
            # norm -1 takes -> min(sum(abs(x), axis=0))
            return np.histogram(r, rBins)[0] * (np.linalg.det(R) / (binVol * len(x1) * len(x2)))
            
        rdf[:,0] += getRDF(x1, x1)
        # rdf[:,1] += getRDF(x1, x2)
        # rdf[:,2] += getRDF(x2, x2)
        nSteps += 1
        
#Plot histgram
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
rdfFile = outFile+".rdf.dat"
plotFile = outFile+".rdf.pdf"
rdf *= (1./nSteps) # normalize
# if rdf[:,0][-1] != 0.:
#     rdf *= (1./rdf[:,0][-1]) # scale to set to 1 far off but doesn't work if 0
np.savetxt(rdfFile, np.hstack((rMid[:,None], rdf)), header='r g11 g12 g22')
plt.plot(rMid, rdf[:,0])
# plt.plot(rMid, rdf[:,1])
plt.xlim(0, rMax)
plt.ylim(0, None)
plt.xlabel('r [A]')
plt.ylabel('g(r)')
plt.legend([legLabel])
# plt.legend(['Na-Na', 'Na-Cl', 'Cl-Cl'])
plt.savefig(plotFile, bbox_inches='tight')
plt.close()

print('DONE')
'''