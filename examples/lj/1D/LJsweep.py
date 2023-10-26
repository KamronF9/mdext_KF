"""Single-site LJ fluid test with roughly the density of water at STP."""
import mdext
import numpy as np
from lammps import PyLammps
from mdext import log
import sys


def main() -> None:
   
    # Current simulation parameters (all in LJ units):
    global T, P
    if sys.argv[1] == -1.:
        P=None
    else:
        P = float(sys.argv[1])  # in LJ epsilon/sigma^2
    T = float(sys.argv[2])  # was 0.7 in LJ epsilon
    seed = int(sys.argv[3])
    U0 = float(sys.argv[4])  # Amplitude of the external potential (in LJ epsilon)
    sigma = 0.5  # Width of the external potential (in LJ sigma)

    # Initialize and run simulation:
    md = mdext.md.MD(
        setup=setup,
        T=T,
        P=P,
        seed=seed,
        potential=mdext.potential.Gaussian(U0, sigma),
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
    md.run(20, "equilibration")  # 10* 50*100 = 50 ps
    md.reset_stats()
    # md.run(40, "collection", f"testPress_{float(sys.argv[1])}.h5")
    md.run(100, "collection", sys.argv[5])


def setup(lmp: PyLammps, seed: int) -> int:
    """Setup initial atomic configuration and interaction potential."""
    
    # Construct simulation box:
    Lz = 30.  # only box dimension that matters
    L = np.array([1., 1., Lz])  # overall box dimensions

    lmp.region(
        f"sim_box block -{L[0]/2} {L[0]/2} -{L[1]/2} {L[1]/2} -{Lz/2} {Lz/2}"
        " units box"
    )
    lmp.create_box("1 sim_box")
    n_bulk = 0.7 # was 0.7   # in LJ 1/sigma

    lmp.log(f'{float(sys.argv[1])}.log')
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
    # lmp.dump(f"write all custom 1000 1D_T{T:.1f}_P{P:.1f}.dump id type x y z")


if __name__ == "__main__":
    main()
