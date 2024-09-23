import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import glob
import matplotlib as mpl
from scipy.ndimage import gaussian_filter1d

sigFilter = 3

def trapz(f: np.ndarray, h: float) -> np.ndarray:
    """Cumulative trapezoidal integral of a function sampled at spacing `h`."""
    return np.concatenate(([0.0], np.cumsum(0.5 * (f[:-1] + f[1:])) * h))

def plot_data() -> None:
    print(glob.glob('*Grouped.h5'))
    Alln = []
    AllE = []
    nLen = []
    for ifname, fname in enumerate(sorted(glob.glob('*Grouped.h5'))):
        # average out the n and E
        print(fname)

        with h5py.File(fname, "r") as f:
            z = np.array(f["z"])
            V = np.array(f["V"])
            lbda = np.array(f["lbda"])
            n = np.array(f["n"])
            E = np.array(f["E"])
            try:
                n_bulk = f.attrs["n_bulk"]
            except:
                n_bulk = None
            # mu = f.attrs["mu"]
            V0 = np.array(f["V0"]) if ("V0" in f) else None
        nLen.append(n.shape[1])
        # print(n.shape) # lambda,grid points
        Alln.append(n.copy())
        # print(len(Alln))
        AllE.append(E.copy())
        # print('AllE',np.array(AllE).shape) #1,lambdas

    # get min len of n
    # clip off uneven data
    cutoff = min(nLen)
    z = z[:cutoff]
    V = V[:cutoff]
    AllnWithCutoff = []
    for i, n in enumerate(Alln):
        print(i)
        # print(np.array(n).shape)
        # print(np.array(n).shape)
        n = np.array(n)[:,:cutoff]
        AllnWithCutoff.append(n)


    # print('np.stack(AllnWithCutoff)',np.stack(AllnWithCutoff).shape) # repl,lambda,points

    # sys.exit(1)
    # print(np.mean(AllnWithCutoff,axis=0).shape) # lbdas, points
    n = np.mean(AllnWithCutoff,axis=0)
    E = np.mean(AllE)

    print(lbda.shape)
    print(lbda)

    # Plot density profiles and potential:
    fig, ax1 = plt.subplots() # was plt.figure(1)
    
    # --- initialize colormap to color by V0
    normalize = mpl.colors.Normalize(vmin=lbda.min(), vmax=lbda.max())
    cmap = mpl.cm.get_cmap("RdBu")
    # --- plot densities
    for lbda_cur, n_cur in zip(lbda, n):
        n_cur = gaussian_filter1d(n_cur,sigFilter)
        ax1.plot(z, n_cur, color=cmap(normalize(lbda_cur)), lw=1) # was plt.
    if n_bulk:
        ax1.axhline(n_bulk, color="k", ls="dotted", lw=1)
    # --- plot potential for comparison
    ax2 = ax1.twinx()
    ax2.plot(z, V, color="k", lw=1, ls="dashed")
    ax1.set_xlabel(r"$z$") 
    ax1.set_ylabel(r"$n(z)$") # was plt.ylabel(r"$n(z)$")
    # ax2.set_ylabel('Potential Shape')
    # ax1.set_xlim(z.min(), z.max())
    ax1.set_xlim([0,None])
    # --- add colorbar
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=normalize)
    sm.set_array([])
    plt.colorbar(sm, label=r"Perturbation strength, $\lambda$", ax=plt.gca())
    plt.savefig(f'_DensityPlots{sigFilter}.pdf')

    if len(lbda) > 1:
        # Compare exact and thermodynamically-integrated energies:
        plt.figure(2)
        # --- exact energies
        if False:
            E -= np.interp(0.0, lbda, E)  # difference from bulk
            plt.plot(lbda, E, label="PE")
        # --- thermodynamic integration
        dz = z[1] - z[0]
        integrand = (n @ V) * dz
        dlbda = lbda[1] - lbda[0]
        E_TI = trapz(integrand, dlbda)
        E_TI -= np.interp(0.0, lbda, E_TI)  # difference from bulk
        plt.plot(lbda, E_TI, "r+", label="TI")
        plt.axhline(0, color="k", lw=1, ls="dotted")
        plt.axvline(0, color="k", lw=1, ls="dotted")
        plt.legend()
        plt.xlim(lbda.min(), lbda.max())
        plt.xlabel(r"Perturbation strength, $V_0$")
        plt.ylabel(r"Free energy change, $\Delta\Phi$")

        # Plot unknown part of potential separately, if applicable:
        if V0 is not None:
            plt.figure()
            for lbda_cur, V0_cur in zip(lbda, V0):
                V_unknown = mu - lbda_cur * V - V0_cur
                plt.plot(z, V_unknown, color=cmap(normalize(lbda_cur)), lw=1)
            plt.xlabel(r"$z$")
            plt.ylabel(r"$V_{\mathrm{unknown}}(z)$")
            plt.xlim(z.min(), z.max())
            # --- add colorbar
            sm = mpl.cm.ScalarMappable(cmap=cmap, norm=normalize)
            sm.set_array([])
            plt.colorbar(sm, label=r"Perturbation strength, $\lambda$", ax=plt.gca())
    plt.savefig('_TIanalysis.pdf')
    plt.show()


def main() -> None:
    # if len(sys.argv) < 2:
    #     print("Usage: python plot.py <data.h5>")
    #     exit(1)
    # data_file = sys.argv[1]
    plot_data()


if __name__ == "__main__":
    main()
