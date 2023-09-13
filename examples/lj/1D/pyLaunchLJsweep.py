import os
from unicodedata import decimal
import numpy as np

# nohup python ../pyLaunchBMH.py > out &
# xx nohup bash ../nacl_bmh.job > out &


# Loop over U in gaussian potential and launch jobs

# s = 1. # sigma width gaussian Ang
# T = 1300. # kelvin
# P = 1. # pressue  can be overridden later in nacl_bmh.py to be none TODO clean up
# p = 2 # atom type to apply the potential to (1-based)
# g = 'planar'

# sweep through potential amplitudes option:
# endRange = 5.0
# stepSize = 0.5

# Ui=-0.4

stepSize = 0.1
    

# for Ui in np.around(np.arange(0,endRange*2 + stepSize ,stepSize), decimals=1)-endRange:  
# direct plug values in option
# for Ui in [-2.5, 2.5]:  
for temperature in np.around(np.arange(1.1,3+stepSize,stepSize),decimals=1): 
    dirName = f'T={temperature}'
    os.system(f'mkdir {dirName}')
    os.chdir(dirName)

    for pressure in np.around(np.arange(0.1,1+stepSize,stepSize),decimals=1):

        # print(f"{Ui:+.1f}")
        # o = f"test-U{Ui:+.1f}.h"
        # log = o[:-3]+"_out"

        # os.system(f"sbatch ../LJsweep.job {Density}") # parallel
        os.system(f"sbatch ../../LJsweep.job {pressure} {temperature}") # parallel

        # serial
        # os.system(f"bash ../LJsweep.job {Density}") 
        # os.system(f'python ../../../parseLammpsLogMdext.py -i log.lammps')
        # os.system(f'mv log.csv {Density}.csv')

        # os.system(f"sbatch ../water_spce.job {Ui} {s} {T} {P} {p} {g} {o}")
        # os.system(f"bash ../nacl_bmh.job {Ui} {s} {T} {P} {p} {g} {o} > {log} &")  
        # break
    
    # os.system(f'bash ../../../parseLammpsLogMdext.sh') # doesn't wait to finish need to run separately
    os.chdir('..')

# os.system('python ../../plotPressVsDensitySingleSaveCSVperT.py')


# THEN run once all are complete and clean up broken runs
# bash ../../../parseLammpsLogMdextManyfolders.sh
# AND
# python ../../../plotPressVsDensitySingleSaveCSVperT.py


   

    
    