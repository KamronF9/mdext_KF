
"""Molten NaCl using Fumi-Tosi potential at 1300 K."""
import mdext
import numpy as np
from lammps import PyLammps
from mdext import MPI, log
from dataclasses import dataclass
import os
import sys


def main(R, startfile, randomSeed, pot) -> None:
    # R = 1.0
    setup = Setup(startfile,pot)
    
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
    pot: str
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

        
        lmp.plugin("load libdeepmd_lmp.so")
        # lmp.pair_style(f'deepmd ../PBED2NaClTrain1AllTraining.pb')
        lmp.pair_style(f'deepmd ../{pot}.pb')
        lmp.pair_coeff("* *")

            # Initial minimize:
        log.info("Minimizing initial structure")
        lmp.minimize("1E-4 1E-6 10000 100000")
    




if __name__ == "__main__":
    
    # main(R)
    # pots = ['D2ClNaPerturbTrain7r11','PBED2NaClTrain1AllTraining']
    # pot = 'D2ClNaPerturbTrain7r11'

    # if len(sys.argv) < 2:
	#     print('Usage: mdext...py <pot file stem wo ext.>')
	#     exit(1)

    # pot = str(sys.argv[1])
    
    if len(sys.argv) < 2:
	    print('Usage: mdext...py R startfile randomSeed pot')
	    exit(1)

    R = float(sys.argv[1])
    startfile = str(sys.argv[2])
    randomSeed = int(sys.argv[3])
    pot = str(sys.argv[4])
    
    main(R,startfile,randomSeed,pot)