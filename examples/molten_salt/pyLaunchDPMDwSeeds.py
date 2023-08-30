import os
from unicodedata import decimal
import numpy as np
import random

# nohup python ../pyLaunchDPMD.py > out &
# nohup bash ../nacl_bmh.job > out &


# Loop over U in gaussian potential and launch jobs

s = 1.0 # sigma width gaussian
T = 1300. # kelvin
P = -1. # -1 is None = NVT, else it takes on the value of P = NPT
p = 2 # atom type to apply the potential to (1-based)
g = 'planar'


# endRange = 5.0
# stepSize = 0.5

# Ui=-0.4
# sweep through potential amplitudes

pots = ['D2ClNaPerturbTrain7r11','PBED2NaClTrain1AllTraining']

for pot in pots:
    # make dir
    os.system(f'mkdir {pot}')
    os.chdir(f'{pot}')
    
    for seedi in range(5):
        
        S = str(random.randint(1,1000))
        Ui = -5.0
        print(f"{Ui:+.1f}_seed{seedi}")
        o = f"test-U{Ui:+.1f}_seed{seedi}.h"
        log = o[:-3]+"_out"

        # orig deepmd
        # os.system(f"bash ../nacl_DPMD_args.job {Ui} {s} {T} {P} {p} {g} {o} {pot}.pb")
        # os.system(f"sbatch ../nacl_DPMD_args.job {Ui} {s} {T} {P} {p} {g} {o} {pot}.pb")
        os.system(f"sbatch ../../nacl_DPMD_args.job {Ui} {s} {T} {P} {p} {g} {o} {pot}.pb {S}")
        # os.system(f"bash ../nacl_bmh.job {Ui} {s} {T} {P} {p} {g} {o}")  #  > {log} &
        # break
   
    os.chdir('..')
    
    