
"""Molten NaCl using SimpleNN trained to PBE+D3 at 1300 K."""
import mdext
import numpy as np
from lammps import PyLammps
from mdext import MPI, log


def main() -> None:

    # Current simulation parameters:
    T = 1300.0  # K
    P = 1.0  # bar
    seed = 12345
    U0 = +3.   # Amplitude of the external potential (kcal/mol)
    sigma = 1. # Width of the external potential (A)

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
        units="metal",
        timestep=0.002,
        Tdamp=0.1,
        Pdamp=0.1,
    )
    md.run(2, "equilibration")
    md.reset_stats()
    md.run(5, "collection", f"test-U{U0:+.0f}.h5")


def setup(lmp: PyLammps, seed: int) -> int:
    """Setup initial atomic configuration and interaction potential."""
    
    # Construct water box:
    L = [20.] * 3  # box dimensions
    file_liquid = "liquid.data"
    is_head = (MPI.COMM_WORLD.rank == 0)
    if is_head:
        mdext.make_liquid.make_liquid(
            pos_min=[-L[0]/2, -L[1]/2, -L[2]/2],
            pos_max=[+L[0]/2, +L[1]/2, +L[2]/2],
            out_file=file_liquid,
            N_bulk=0.015,
            masses=[22.99, 35.45],
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
    lmp.pair_coeff("1 1 0.2637 0.317 2.340 1.048553 -0.49935") # Na-Na
    lmp.pair_coeff("2 2 0.158221 0.327 3.170 75.0544 -150.7325")  # Cl-Cl
    lmp.pair_coeff("1 2 0.21096 0.317 2.755 6.99055303 -8.6757")  # Na-Cl
    lmp.kspace_style("ewald 1e-5")

    # Initial minimize:
    log.info("Minimizing initial structure")
    lmp.minimize("1E-4 1E-6 10000 100000")


if __name__ == "__main__":
    main()
