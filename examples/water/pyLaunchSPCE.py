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
endRange = 2.5
stepSize = 0.5

# Ui=-0.4


    

for Ui in np.around(np.arange(0,endRange*2 + stepSize ,stepSize), decimals=1)-endRange:  
# direct plug values in option
# for Ui in [-2.5, 2.5]:  
    print(f"{Ui:+.1f}")
    # o = f"test-U{Ui:+.1f}.h"
    # log = o[:-3]+"_out"

    os.system(f"sbatch ../spce.job {Ui}")
    # os.system(f"sbatch ../water_spce.job {Ui} {s} {T} {P} {p} {g} {o}")
    # os.system(f"bash ../nacl_bmh.job {Ui} {s} {T} {P} {p} {g} {o} > {log} &")  
    # break
   

    
    