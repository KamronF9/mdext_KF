import os
from unicodedata import decimal
import numpy as np
import random

# nohup python ../pyLaunchBMH.py > out &
# xx nohup bash ../nacl_bmh.job > out &


# Loop over U in gaussian potential and launch jobs
Ts = [0.7, 1.0, 1.6]
Ps = [0.7, 1.0, 1.3]

for i in range(len(Ts)):
    # i = 0
    
    # s = 1.0 # sigma width gaussian
    T = Ts[i] # LJ epsilon units
    P = Ps[i] # -1 is None = NVT, else it takes on the value of P = NPT
    # p = 1 # atom type to apply the potential to (1-based)
    # g = 'planar'
    os.system('cp ../LJsweep.job ../LJsweep.py ../pyLaunchLJExtPotSweep.py .')
    
    newDir = f'P{P}_T{P}'
    os.system(f'mkdir {newDir}')
    os.chdir(newDir)
    os.system('rm out.o*')



    # sweep through potential amplitudes option:
    # endRange = 15.0
    # stepSize = 0.5

    # Ui=-0.4

    # stepSize = 0.5



    # direct plug values in option
    # for Ui in [-2.5, 2.5]:  
    # for temperature in np.around(np.arange(1.1,3+stepSize,stepSize),decimals=1): 
    # for temperature in np.around(np.arange(0.1,1.5+stepSize,stepSize),decimals=1): 
    #     dirName = f'T={temperature}'
    #     os.system(f'mkdir {dirName}')
    #     os.chdir(dirName)
    # for pressure in np.around(np.arange(0.1,1+stepSize,stepSize),decimals=1):

    # seeds = np.arange(1,2,1) # seed for np in the launches, doesn't include end
    # seeds = np.arange(1,3,1) # seed for np in the launches, doesn't include end
    seeds = [0]

    for seed in seeds:
        print('seed: ',seed)
        os.system(f'mkdir seed{seed:04}')
        os.chdir(f'seed{seed:04}')
        os.system('rm out.o*')
        
        os.system(f"sbatch ../../../LJsweep.job {P} {T} {seed} pos") # parallel
        os.system(f"sbatch ../../../LJsweep.job {P} {T} {seed} neg") # parallel
        os.chdir('..')


    os.chdir('..')
