import os
from unicodedata import decimal
import numpy as np
import random

# seeds = np.arange(1,2,1) # seed for np in the launches, doesn't include end
# seeds = np.arange(1,3,1) # seed for np in the launches, doesn't include end
# seeds = [0,1,2]

# launch in folder for all data
# copy files in
os.system('cp ../05* .')

replicas = np.arange(1,5,1) # 4 replicas

# seeds = np.arange(0,16,1) # seed for np in the launches, doesn't include end
# seeds = [1]
totalSeedsToRun = 48 # if starting from 0 in arange
seeds = np.arange(0,totalSeedsToRun,1) # seed for np in the launches, doesn't include end

for seed in seeds:
    print('seed: ',seed)
    seedName = f'seed{seed:04}'
    os.system(f'mkdir {seedName}')
    os.chdir(seedName)
    
    for replica in replicas:
        replName = f'replica{replica:04}'
        os.system(f'mkdir {replName}')
        os.chdir(replName)

        os.system('rm out.o*')
        os.system(f'sbatch ../../05spceMono.job {seed}') # parallel
        # os.system(f"sbatch ../../../LJsweep.job {P} {T} {seed} neg") # parallel
        os.chdir('..')


    os.chdir('..')


os.chdir('..')
