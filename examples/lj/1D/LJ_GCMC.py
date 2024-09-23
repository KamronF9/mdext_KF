"""Single-site LJ fluid test with roughly the density of water at STP."""
import mdext
import numpy as np
from lammps import PyLammps
from mdext import log
import sys
import scipy
from dataclasses import dataclass
import os

def main(P, T, U0, startfile, npSeed) -> None:
    mu = 500.
    displace = 0.1 #0 1.
    overlap_cutoff = 0.2 # LJ distance units - sigma
    
    # Ptarget = 0.7 # LJ units epsilon/sigma^2
    sigma = 1.0  # Width of the external potential (in LJ sigma)

    np.random.seed(npSeed)
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
        seed=12345,
        potential=mdext.potential.Gaussian(U0, sigma, coeffs),
        geometry_type=mdext.geometry.Planar,
        n_atom_types=1,
        potential_type=1,
        pe_collect_interval=100,
        dimension=1,
        units="lj",
        timestep=0.01,
        Tdamp=0.5,
        Pdamp=1.0,
    )
    # GCMC
    # fix ID group-ID gcmc N X M type seed T mu displace keyword values
    # region sim_box pressure {Ptarget} overlap_cutoff {overlap_cutoff}
    md.lmp.fix(f'gcmc1 all gcmc 100 3 3 1 1234 {T} {mu} {displace} region atom_box full_energy')
    md.run(20, "equilibration")  # 10* 50*100 = 50 ps
    # md.reset_stats()
    # md.run(40, "collection", f"testPress_{float(sys.argv[1])}.h5")
    md.run(100, "collection", 'test.h5') # f"data-U{U0:+.2f}.h5")
    md.lmp.write_data(f'U{U0:+.2f}.step.data nocoeff')


@dataclass
class Setup:
    startfile: str

    def __call__(self, lmp: PyLammps, seed: int) -> int:
        """Setup initial atomic configuration and interaction potential."""
        
        file_liquid = self.startfile


        # Construct simulation box:
        Lz = 30.  # only box dimension that matters
        L = np.array([1., 1., Lz])  # overall box dimensions

        lmp.region(
            f"sim_box block -{L[0]/2} {L[0]/2} -{L[1]/2} {L[1]/2} -{Lz/2} {Lz/2}"
            " units box"
        )
        lmp.create_box("1 sim_box")

        # lmp.region(
        #     f"mc_box block -{L[0]/100000} {L[0]/100000} -{L[1]/100000} {L[1]/100000} -{Lz/2} {Lz/2}"
        #     " units box"
        # )
        n_bulk = 0.7 # was 0.7   # in LJ 1/sigma

        # lmp.log(f'{float(sys.argv[1])}.log')
        n_atoms = int(np.round(n_bulk * Lz))
        lmp.region(f"atom_box block -0.0 0.0 -0.0 0.0 -{Lz/2} {Lz/2} units box")
        lmp.create_atoms(f"1 random {n_atoms} {seed} atom_box")
        lmp.mass("1 1.")

        # Interaction potential:
        lmp.pair_style("lj/cut 10")
        lmp.pair_coeff("1 1 1.0 1.0")
        # lmp.pair_modify("tail yes")

        # Initial minimize:
        log.info("Minimizing initial structure")
        lmp.minimize("1E-4 1E-6 10000 100000")

        # Dump output file - conflicts with the reset stats function above somehow
        lmp.dump(f"write all custom 100 1D_GCMC.dump id type x y z")
        # lmp.dump(f"write all custom 1000 1D_T{T:.1f}_P{P:.1f}.dump id type x y z")


if __name__ == "__main__":
    
    # if (len(sys.argv)-1) < 4:
    #     print('usage python ...py <P> <T> <seed> <sign of potentials>')
    #     sys.exit(1)

    # npSeed = 1 # np seed for generating the potential shape 
    # each seed would be run in a separate folder to prevent conflict and run in parallel
    
    # Current simulation parameters (all in LJ units):
    os.system('rm log.lammps 1D_GCMC.dump')
    P=None # NVT
    # if sys.argv[1] == -1.:
    #     P=None
    # else:
    #     P = float(sys.argv[1])  # in LJ epsilon/sigma^2
    T = 0.7  # was 0.7 in LJ epsilon
    npSeed = 2
    
    posOrNegLbdas = 'pos' # sys.argv[4]  # pos or neg
    # sweep through potential amplitudes option:

    endRange = 15.0 # Amplitude of the external potential (in LJ epsilon)
    stepSize = 0.5 
    
    # Ui in eV
    UisPos = np.around(np.arange(0,endRange + stepSize ,stepSize), decimals=2) # positive only
    UisNeg = np.around(np.arange(0,-endRange - stepSize ,-stepSize), decimals=2) # Negative
    # print('Uis: ', UisPos)
    # print('Uis: ', UisNeg)
    # Ui=-0.4
    # sys.exit(1)
    # for Uis in [UisPos,UisNeg]:
    if posOrNegLbdas=='pos':
        Uis = UisPos
    elif posOrNegLbdas=='neg':
        Uis = UisNeg
    else:
        print('invalid argument sign of potential')
        sys.exit(1)

    print('Uis to run: ', Uis)

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
        main(P,T,Ui,startfile,npSeed)

        break
        # while not (os.path.exists(pwd+'/'+currdatafile)):
            # time.sleep(1)   
        # if i==1: break  # HACK

    print('Done!')