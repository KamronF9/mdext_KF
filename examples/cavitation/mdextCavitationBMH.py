
"""Molten NaCl using Fumi-Tosi potential at 1300 K."""
import mdext
import numpy as np
from lammps import PyLammps
from mdext import MPI, log
from dataclasses import dataclass


def main(R, startfile, randomSeed) -> None:
    # R = 1.0
    setup = Setup(startfile)
    
    # Initialize and run simulation:
    md = mdext.md.MD(
        setup=setup,
        T=1100.0,
        P=1.0,
        seed=randomSeed,
        potential=mdext.potential.Exponential(A=0.0211, rho=0.25, sigma=R),
        geometry_type=mdext.geometry.Spherical,
        n_atom_types=2,
        potential_type=0,
        pe_collect_interval=50,
        units="metal",
        timestep=0.002,
        Tdamp=0.1,
        Pdamp=0.1,
    )

    md.run(2, "equilibration")
    md.reset_stats()
    md.run(5, "collection", f"{R}_cavitation.h5")
    md.lmp.write_data(str(R) + '.cavity.data nocoeff')

@dataclass
class Setup:
    startfile: str
    # R: float = 0.0
    
    # replaces 
    # def __init__(self, startfile):
    #     self.startfile = startfile
    
    def __call__(self, lmp: PyLammps, seed: int) -> int:
        """Setup initial atomic configuration and interaction potential."""

        # Construct water box:
        L = [20.] * 3  # box dimensions

        file_liquid = self.startfile

        is_head = (MPI.COMM_WORLD.rank == 0)
        if is_head and file_liquid == 'liquid.data':
            mdext.make_liquid.make_liquid(
                pos_min=[-L[0]/2, -L[1]/2, -L[2]/2],
                pos_max=[+L[0]/2, +L[1]/2, +L[2]/2],
                out_file=file_liquid,
                N_bulk=0.015,
                masses=[35.45, 22.99],
                radii=[1.0, 1.0],
                atom_types=[1, 2],
                atom_pos=[[0., 0., 1.3], [0., 0., -1.3]],
                bond_types=np.zeros((0,), dtype=int),
                bond_indices=np.zeros((0, 2), dtype=int),
                angle_types=np.zeros((0,), dtype=int),
                angle_indices=np.zeros((0, 3), dtype=int),
            )

        lmp.atom_style("full")
        lmp.read_data(file_liquid)
        # Interaction potential (Fumi-Tosi w/ Ewald summation):
        lmp.pair_style("born/coul/long 9.0")
        lmp.pair_coeff("2 2 0.2637 0.317 2.340 1.048553 -0.49935") # Na-Na
        lmp.pair_coeff("1 1 0.158221 0.327 3.170 75.0544 -150.7325")  # Cl-Cl
        lmp.pair_coeff("1 2 0.21096 0.317 2.755 6.99055303 -8.6757")  # Na-Cl
        lmp.set("type 2 charge +1")
        lmp.set("type 1 charge -1")
        lmp.kspace_style("pppm 1e-5")

            # Initial minimize:
        log.info("Minimizing initial structure")
        lmp.minimize("1E-4 1E-6 10000 100000")
    




if __name__ == "__main__":
    
    # main(R)
    
    # run cavitation sims
    # batch them 
    #Rs = np.arange(0.2,9.1,0.1) ORIG
    #Rs = np.arange(0.0,0.3,0.1)
    Rs = np.arange(0.0,9.1,0.1)
    Rs = list(Rs)
    Rs = [round(R,1) for R in Rs]
    #print(Rs)
    #exit()
    for i, R in enumerate(Rs):
        #print(R)
        #break
        R = round(R,1)
        print(f"launching cavity {R}")
        #os.system("lmp_0921 < cavity.in -v T 1100 -v R %s"%(str(R)))
        if R == 0.0:  # initial run
            startfile = 'liquid.data'
    #    if R == 8.3:  # continue initial run
    #        startfile = '8.2.cavity.data'
        else:  # normal sequence within a run
            prev_R = Rs[Rs.index(R)-1]
            startfile = str(prev_R)+'.cavity.data'  
        randomSeed = np.random.randint(0,1000)
        # switch if dpmd or cavity
        # pot='PBED2NaClTrain1AllTraining'
        # os.system("nohup bash cavityDPMD.job %s %s %s %s > LOG &"%(str(R),startfile,str(randomSeed),pot))
        # os.system("nohup bash cavity.job %s %s %s > LOG &"%(str(R),startfile,str(randomSeed)))
        main(R,startfile,randomSeed)
        # break
        # currdatafile = str(R) + '.cavity.data'
        # pwd = os.getcwd()
        # while not (os.path.exists(pwd+'/'+currdatafile)):
            # time.sleep(1)   
        # if i==1: break  # HACK
    print('Done!')