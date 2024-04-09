from ase import units
from ase.md.langevin import Langevin
from ase.io import read, write
from ase.io.lammpsdata import read_lammps_data
import numpy as np
import time
import sys
import os
from ase import Atom, Atoms
from ase.build import bulk

from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from ase.calculators.lammpsrun import LAMMPS
from deepmd.calculator import DP
from deepmd.infer import DeepPot
from ase.calculators.mixing import SumCalculator

modelTypes = ['deepwater','NNP', 'NNPext', 'DFT']
models = ['/global/homes/k/kamron/Scratch/NNmd/water/ice/Ice-XI-test/H2O-Phase-Diagram-model_compressed.pb','/global/homes/k/kamron/Scratch/AIMDtoNN/Water/Train1Reg/seed0baseline/H1O2WaterTrain1r6_168k.pb','/global/homes/k/kamron/Scratch/AIMDtoNN/Water/Train4PerturbW_O_LVint/BaselineR0_works/H1O2WaterPerturbTrain4wO_LVintSeed0.pb','']



frames = read('water_t300_p1.xyz', ':') 
# frames = read('T300_surf.xyz', ':') 

# remove first half of frames
halfPoint = int(len(frames)/2)
frames = frames[halfPoint:]

energies = np.zeros([len(modelTypes),len(frames)])


for imodel,modelType in enumerate(modelTypes):
    

    model_path = models[imodel]
    if modelType != 'DFT':
        dp=DeepPot(model_path)


    # deepwater=False  #

    # datain = '/global/homes/k/kamron/Scratch/NNmd/water/ice/Ice-XI-test/Bulk/iceBulk.data'
    # ats = read_lammps_data(datain)
    # print(frames.get_chemical_symbols())
    # atoms.arrays["bonds"]

    for iframe, ats in enumerate(frames):
        print(iframe)
        # print(ats.get_chemical_symbols())
        # atype = (np.array(ats.todict()['type'].tolist())-1).tolist()
        atype = ats.get_atomic_numbers()

        # mapping is H->1 O->8 in xyz read in but need to cast in right mode for deepmd
        if modelType=='deepwater':
            # but need to map to H 1 O 0 for deepwater
            # tempType = np.where(ats.get_atomic_numbers()==1,2,ats.get_atomic_numbers()) 
            # atype = np.where(tempType==8,2,ats.get_atomic_numbers()) 
            d = {8:0,1:1}
            # d = {8:1,1:0}
        else:
            # need to map to H 0 O 1 for NNP NNPext
            # atype = np.where(ats.get_atomic_numbers()==8,2,ats.get_atomic_numbers()).tolist()
            # 1 is fine as H
            d = {8:1,1:0}
        atype_translated = [d[i] for i in atype]
        
        


        # atype = [1,0,1]
        # e, f, v = dp.eval(coord, cell, atype)
        
        # atype = (np.array(ats.todict()['type'].tolist())-1).tolist()
        # e, f, v = dp.eval(ats.get_positions().flatten(),ats.get_cell().flatten(),atype)
        if modelType == 'DFT':
            e = ats.get_total_energy()
            # print(e)
            # energies[imodel,iframe]=e
        else:
            e, f, v = dp.eval(ats.get_positions().flatten(),ats.get_cell().flatten(),atype_translated)
            # print(e[0][0])
            e=e[0][0]

        energies[imodel,iframe]=e
        
        # if iframe==1: break
        # break

print(energies)
print(np.mean(energies,axis=1))