#!/usr/bin/env python
# RDF only
import numpy as np
import sys

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
