"""SPC/E water test at STP."""
import mdext
import numpy as np
from lammps import PyLammps
from mdext import MPI, log
import sys
import scipy
from dataclasses import dataclass
import random

def main(U0ev, startfile, npSeed) -> None:

    # Current simulation parameters:
    T = 300.  # K
    P = 1  # atm
    # P = None  # NVT
    seed = random.randint(1,1000) # mdext lammps seed
    # seed = 345 # mdext lammps seed last ran with 
    # global U0ev
    # U0ev = (float(sys.argv[1]))  # Amplitude of the external potential (eV)
    U0 = U0ev * 23.06   # Convert amplitude of the external potential to kcal/mol since units water classical 
    
    # sigma = 1.0 # Width of the external potential (A)
    # np.random.seed(npSeed)
    # coeffs = np.random.randn(3) # gauss and poly
    # powers = np.arange(len(coeffs))
    # power_pair_sums = powers[:, None] + powers[None, :]
    # norm_fac = np.sqrt(
    #     np.sqrt(np.pi)
    #     / (coeffs @ scipy.special.gamma(power_pair_sums + 0.5) @ coeffs)
    # )
    # coeffs *= norm_fac


    # U0 = 1.  # lambda equivalent to scale all
    np.random.seed(npSeed)
    sigma = np.random.uniform(0.5, 3)  # ang
    B = np.random.randn() 
    if B > 0.5:
        A = -1.
    elif B < 0.:
        A = 1.
    else:
        A = np.random.choice([1,-1]) 

    coeffs = A*np.array([1, B])

    
    setup = Setup(startfile)

    # Initialize and run simulation:
    md = mdext.md.MD(
        setup=setup,
        T=T,
        P=P,
        seed=seed,
        potential=mdext.potential.Gaussian(U0, sigma, coeffs),
        geometry_type=mdext.geometry.Planar,
        n_atom_types=2,
        potential_type=2,
        pe_collect_interval=100,
        # Pdamp=1000.,
        # timestep=1,
    )
    # save dump then remove reset stats
    # md.lmp.dump(f"write all custom 1000 pylammps{U0ev:+.2f}.dump id type element x y z")
    # md.lmp.dump_modify("write element H O")
    md.run(5, "equilibration") # was 5 . 5000 dt per run #
    md.reset_stats() # was commented to allow dump?
    md.run(20, "collection", f"data-U{U0ev:+.2f}.h5") # was 10
    md.lmp.write_data(f'U{U0ev:+.2f}.step.data nocoeff')


@dataclass
class Setup:
    startfile: str

    def __call__(self, lmp: PyLammps, seed: int) -> int:
        """Setup initial atomic configuration and interaction potential."""
        
        file_liquid = self.startfile

        # Construct water box:
        L = np.array([30.,30.,40.])  # overall box dimensions
        # L = np.array([26.44, 24.67, 23.60])  # overall box dimensions
        # L = np.array([26.44/2, 24.67/2, 23.60/2])  # overall box dimensions
        # file_liquid = "liquid.data"
        is_head = (MPI.COMM_WORLD.rank == 0)
        if is_head and file_liquid == 'liquid.data':
            mdext.make_liquid.make_water(
                pos_min=[-L[0]/2, -L[1]/2, -L[2]/2],
                pos_max=[+L[0]/2, +L[1]/2, +L[2]/2],
                out_file=file_liquid,
            )
        lmp.atom_style("full")
        lmp.read_data(file_liquid)

        # Interaction potential (SPC/E, long-range):
        lmp.pair_style("lj/cut/coul/cut 10")
        lmp.bond_style("harmonic")
        lmp.angle_style("harmonic")
        # lmp.kspace_style("pppm/disp 1e-5")
        # lmp.kspace_modify("mesh/disp 5 5 5 gewald/disp 0.24")
        lmp.set("type 1 charge  0.4238")
        lmp.set("type 2 charge -0.8476")
        lmp.pair_coeff("1 *2 0.000 0.000")  # No LJ for H
        lmp.pair_coeff("2 2 0.1553 3.166")  # O-O
        lmp.bond_coeff("1 1000 1.0")  # H-O
        lmp.angle_coeff("1 100 109.47")  # H-O-H

        # Initial minimize:
        log.info("Minimizing initial structure")
        lmp.neigh_modify("exclude molecule/intra all")
        lmp.minimize("1E-4 1E-6 10000 100000")
        
        # Rigid molecule constraints for dynamics:
        
        lmp.neigh_modify("exclude none")
        # lmp.dump("write all custom 100 20_2.lammpstrj id type x y z vx vy vz")
        lmp.fix("BondConstraints all shake 0.001 20 0 b 1 a 1")



if __name__ == "__main__":
    
    npSeed = int(sys.argv[1]) # np seed for generating the potential shape
    # npSeed = 1 # np seed for generating the potential shape 
    # each seed would be run in a separate folder to prevent conflict and run in parallel
    # posOrNegLbdas = sys.argv[1]  # pos or neg
    # sweep through potential amplitudes option:
    endRange = 0.2 #1.0 # eV
    stepSize = 0.02 #0.01  # start around kbT
    
    # Ui in eV
    UisPos = np.around(np.arange(0,endRange + stepSize ,stepSize), decimals=2) # positive only
    # UisNeg = np.around(np.arange(0,-endRange - stepSize ,-stepSize), decimals=2) # Negative
    # print('Uis: ', UisPos)
    # print('Uis: ', UisNeg)
    # Ui=-0.4
    # sys.exit(1)
    # for Uis in [UisPos,UisNeg]:
    # if posOrNegLbdas=='pos':
    Uis = UisPos
    # elif posOrNegLbdas=='neg':
    #     Uis = UisNeg
    # else:
    #     print('invalid argument 1')
    #     sys.exit(1)

    print('Uis to run: ', Uis)

    for i, Ui in enumerate(Uis):  
    # direct plug values in option
    # for Ui in [-2.5, 2.5]:  
        # print(f"{Ui:+.1f}")
        
        
        #os.system("lmp_0921 < cavity.in -v T 1100 -v R %s"%(str(R)))
        if Ui == 0.0:  # initial run
            # use expanded cell file
            # startfile = 'U+0.00.step.data' # skip now that we already ran equilibrium 
            startfile = 'liquid.data' # hack 
            #    if R == 8.3:  # continue initial run
            #        startfile = '8.2.cavity.data'
        else:  # normal sequence within a run
            # prev_Ui = UisPos[UisPos.index(Ui)-1]
            prev_Ui = Uis[i-1]
            startfile = f'U{prev_Ui:+.2f}.step.data'  
        # randomSeed = np.random.randint(0,1000)

        print(f"launching seed {npSeed}, Ui {Ui:+.2f}")
        main(Ui,startfile,npSeed)

        # break
        # while not (os.path.exists(pwd+'/'+currdatafile)):
            # time.sleep(1)   
        # if i==1: break  # HACK

    print('Done!')


# endRange = 0.1 # eV
# stepSize = 0.01  # start around kbT
# UisPos=np.around(np.arange(0,endRange + stepSize ,stepSize), decimals=2).tolist()
# UisPos[UisPos.index(0.02)-1]
# for i, Ui in enumerate(UisPos): 
#     print('i ',UisPos[i])
#     if i!=0:
#         print('i-1 ',UisPos[i-1])