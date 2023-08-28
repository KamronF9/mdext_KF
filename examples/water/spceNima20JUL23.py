"""SPC/E water test at STP."""
import mdext
import numpy as np
from lammps import PyLammps
from mdext import MPI, log
import sys

def main() -> None:

    # Current simulation parameters:
    T = 298.0  # K
    P = None  # atm
    seed = 345s
    U0 = (float(ys.argv[1]))*23.06   # Amplitude of the external potential (kcal/mol)
    sigma = 1.0 # Width of the external potential (A)
    global jj
    jj=float(sys.argv[1]) # eV
    # Initialize and run simulation:
    md = mdext.md.MD(
        setup=setup,
        T=T,
        P=P,
        seed=seed,
        potential=mdext.potential.Gaussian(U0, sigma),
        geometry_type=mdext.geometry.Planar,
        n_atom_types=2,
        potential_type=2,
        timestep=1,
    )
    md.run(1, "equilibration")
    md.reset_stats() # was commented to allow dump?
    md.run(2, "collection", f"water-lj-pl-U{jj:+.1f}-{sigma:+.2f}.h5")


def setup(lmp: PyLammps, seed: int) -> int:
    """Setup initial atomic configuration and interaction potential."""
    
    # Construct water box:
    # L = np.array([26.44, 24.67, 23.60])  # overall box dimensions
    L = np.array([26.44/2, 24.67/2, 23.60/2])  # overall box dimensions
    file_liquid = "liquid.data20"
    is_head = (MPI.COMM_WORLD.rank == 0)
    if is_head:
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
    lmp.dump(f"write all custom 1000 pylammps{jj:+.1f}.dump id type element x y z")
    lmp.dump_modify("write element H O")


if __name__ == "__main__":
    main()