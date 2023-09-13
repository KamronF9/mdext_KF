import os
from unicodedata import decimal
import numpy as np

# nohup python ../pyLaunchDPMD.py > out &
# nohup bash ../nacl_bmh.job > out &


# Loop over U in gaussian potential and launch jobs

s = 1.0 # sigma width gaussian
T = 1300. # kelvin
P = -1. # -1 is None = NVT, else it takes on the value of P = NPT
p = 2 # atom type to apply the potential to (1-based)
g = 'planar'
S = 1234 # seed


# endRange = 5.0
# stepSize = 0.5

# Ui=-0.4
# sweep through potential amplitudes

# pots = ['D2ClNaPerturbTrain7r11','PBED2NaClTrain1AllTraining']
# pots = ['PBED2NaClTrain8Pert_noHighPot']
pots = []
for i in range(4):
    pots.append('PBED2NaClTrain8Pert_noHighPotseed'+str(i))
    pots.append('PBED2NaClTrain1AllTrainingseed'+str(i))

for pot in pots:
    # make dir
    os.system(f'mkdir {pot}')
    os.chdir(f'{pot}')
    
    for Ui in [-5.0, 5.0]:      
    # for Ui in np.around(np.arange(0,endRange*2 + stepSize ,stepSize), decimals=1)-endRange:  
        print(f"{Ui:+.1f}")
        o = f"test-U{Ui:+.1f}.h"
        log = o[:-3]+"_out"

        # orig deepmd
        # os.system(f"bash ../nacl_DPMD_args.job {Ui} {s} {T} {P} {p} {g} {o} {pot}.pb")
        # os.system(f"sbatch ../nacl_DPMD_args.job {Ui} {s} {T} {P} {p} {g} {o} {pot}.pb")
        os.system(f"sbatch ../../nacl_DPMD_args.job {Ui} {s} {T} {P} {p} {g} {o} {pot}.pb {S}")
        # os.system(f"bash ../nacl_bmh.job {Ui} {s} {T} {P} {p} {g} {o}")  #  > {log} &
        # break
   
    os.chdir('..')
    
    