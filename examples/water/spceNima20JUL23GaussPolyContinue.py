"""SPC/E water test at STP."""
import mdext
import numpy as np
from lammps import PyLammps
from mdext import MPI, log
import sys
import scipy
from dataclasses import dataclass

def main(U0ev, startfile, npSeed) -> None:

    # Current simulation parameters:
    T = 298.0  # K
    P = 1  # atm
    # P = None  # atm
    seed = 345
    # global U0ev
    # U0ev = (float(sys.argv[1]))  # Amplitude of the external potential (eV)
    U0 = U0ev * 23.06   # Amplitude of the external potential (kcal/mol)
    sigma = 1.0 # Width of the external potential (A)
    np.random.seed(0)
    coeffs = np.random.randn(3) # gauss and poly
    powers = np.arange(len(coeffs))
    power_pair_sums = powers[:, None] + powers[None, :]
    norm_fac = np.sqrt(
        np.sqrt(np.pi)
        / (coeffs @ scipy.special.gamma(power_pair_sums + 0.5) @ coeffs)
    )
    coeffs *= norm_fac
    
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
        timestep=1,
    )
    md.run(5, "equilibration") # was 5 . 5000 dt per run #
    md.reset_stats() # was commented to allow dump?
    md.run(10, "collection", f"data-U{U0ev:+.2f}.h5") # was 10
    md.lmp.write_data(f'U{U0ev:+.2f}.step.data nocoeff')

@dataclass
class Setup:
    startfile: str

    def __call__(self, lmp: PyLammps, seed: int) -> int:
        """Setup initial atomic configuration and interaction potential."""
        
        file_liquid = self.startfile

        # Construct water box:
        L = np.array([26.44, 24.67, 23.60])  # overall box dimensions
        # L = np.array([26.44/2, 24.67/2, 23.60/2])  # overall box dimensions
        file_liquid = "liquid.data"
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
        lmp.pair_style("lj/long/coul/long long long 10")
        lmp.bond_style("harmonic")
        lmp.angle_style("harmonic")
        lmp.kspace_style("pppm/disp 1e-5")
        lmp.kspace_modify("mesh/disp 5 5 5 gewald/disp 0.24")
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
        # # lmp.fix("BondConstraints all shake 0.001 20 0 b 1 a 1")
        # save dump then remove reset stats
        if False:
            lmp.dump(f"write all custom 1000 pylammps{U0ev:+.1f}.dump id type element x y z")
            lmp.dump_modify("write element H O")


if __name__ == "__main__":
    
    # npSeed = int(sys.argv[1]) # np seed for generating the potential shape
    npSeed = 1 # np seed for generating the potential shape 
    # each seed would be run in a separate folder to prevent conflict and run in parallel

    # sweep through potential amplitudes option:
    endRange = 1.0 # eV
    stepSize = 0.01  # start around kbT

    # Ui in eV
    UisPos = np.around(np.arange(0,endRange + stepSize ,stepSize), decimals=2) # positive only
    UisNeg = np.around(np.arange(0,-endRange - stepSize ,-stepSize), decimals=2) # Negative
    print('Uis: ', UisPos)
    print('Uis: ', UisNeg)
    # Ui=-0.4
    # sys.exit(1)
    for Uis in [UisPos,UisNeg]:

        for i, Ui in enumerate(Uis):  
        # direct plug values in option
        # for Ui in [-2.5, 2.5]:  
            # print(f"{Ui:+.1f}")
            
            print(f"launching seed {npSeed}, Ui {Ui:+.2f}")
            #os.system("lmp_0921 < cavity.in -v T 1100 -v R %s"%(str(R)))
            if Ui == 0.0:  # initial run
                startfile = 'liquid.data'
        #    if R == 8.3:  # continue initial run
        #        startfile = '8.2.cavity.data'
            else:  # normal sequence within a run
                # prev_Ui = UisPos[UisPos.index(Ui)-1]
                prev_Ui = Uis[i-1]
                startfile = f'U{prev_Ui:+.2f}.step.data'  
            # randomSeed = np.random.randint(0,1000)
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