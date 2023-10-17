import os
from unicodedata import decimal
import numpy as np
import random

# nohup python ../pyLaunchBMH.py > out &
# xx nohup bash ../nacl_bmh.job > out &


# Loop over U in gaussian potential and launch jobs
Ts = [1.0, 1.3, 1.6]
Ps = [0.7, 0.7, 1.0]
i = 0

# s = 1.0 # sigma width gaussian
T = Ts[i] # LJ epsilon units
P = Ps[i] # -1 is None = NVT, else it takes on the value of P = NPT
# p = 1 # atom type to apply the potential to (1-based)
# g = 'planar'

# sweep through potential amplitudes option:
endRange = 15.0
stepSize = 0.5

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

# for Ui in np.around(np.arange(0,endRange*2 + stepSize ,stepSize), decimals=1)-endRange:  
for Ui in [12.5]:  
    print(f"{Ui:+.1f}")
    h5name = f"test-U{Ui:+.1f}.h"
    seed = str(random.randint(1,1000))
    # log = o[:-3]+"_out"

    # os.system(f"sbatch ../LJsweep.job {Density}") # parallel
    os.system(f"sbatch ../LJsweep.job {P} {T} {seed} {Ui} {h5name}") # parallel

 
    # break
