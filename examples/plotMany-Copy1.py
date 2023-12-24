import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import os
import matplotlib as mpl

endRange = 5.0
stepSize = 0.5

particle = 1

plt.figure(1)
# --- initialize colormap to color by V0
normalize = mpl.colors.Normalize(vmin=-endRange, vmax=endRange)
cmap = mpl.cm.get_cmap("RdBu")

for Ui in np.around(np.arange(0,endRange*2 + stepSize ,stepSize), decimals=1)-endRange:  
    print(f"{Ui:.1f}")
    
    with h5py.File(f"test-U-{Ui:.1f}.h5", "r") as fp:
        r = np.array(fp["r"])
        n = np.array(fp["n"])
        V = np.array(fp["V"])
    plt.plot(r,n[:,particle]/(1.),color=cmap(normalize(Ui)), lw=1)
    plt.xlabel("z [$\AA$]",fontsize=11,fontname="Times New Roman")
    plt.ylabel("$n_{H}(z)$",fontsize=11,fontname="Times New Roman")
    # plt.ylabel("$n_{H}(z)/n_{bulk}$",fontsize=11,fontname="Times New Roman")
    plt.plot(r, V, color="k", lw=1, ls="dashed")  # plot largest shaped V
    # plt.legend()
    plt.legend( prop={
            'family': 'Times New Roman', "size": 7, 'stretch': 'normal'})    # plt.title('Classical, -0.4 eV',fontsize=12,fontname="Times New Roman")
    plt.xticks(fontsize=11,fontname="Times New Roman")
    plt.yticks(fontsize=11,fontname="Times New Roman")


plt.xlim([r.min(),r.max()])
# plt.ylim([0,4])

# --- add colorbar
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=normalize)
sm.set_array([])
plt.colorbar(sm, label=r"Perturbation strength, $\lambda$")
figure = plt.gcf()
figure.set_size_inches(3.35, 2.2)
plt.savefig('plotMany.pdf', dpi=600, bbox_inches='tight')



def plot_data(data_file: str) -> None:
    f = h5py.File(data_file, "r")
    z = np.array(f["z"])
    V = np.array(f["V"])
    lbda = np.array(f["lbda"])
    n = np.array(f["n"])
    E = np.array(f["E"])
    n_bulk = f.attrs["n_bulk"]

    # Plot density profiles and potential:
    plt.figure(1)
    # --- initialize colormap to color by V0
    normalize = mpl.colors.Normalize(vmin=lbda.min(), vmax=lbda.max())
    cmap = mpl.cm.get_cmap("RdBu")
    # --- plot densities
    for lbda_cur, n_cur in zip(lbda, n):
        plt.plot(z, n_cur, color=cmap(normalize(lbda_cur)), lw=1)
    plt.axhline(n_bulk, color="k", ls="dotted", lw=1)
    # --- plot potential for comparison
    plt.plot(z, V, color="k", lw=1, ls="dashed")
    plt.xlabel(r"$z$")
    plt.ylabel(r"$n(z)$")
    plt.xlim(z.min(), z.max())
    # --- add colorbar
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=normalize)
    sm.set_array([])
    plt.colorbar(sm, label=r"Perturbation strength, $\lambda$")

    if len(lbda) > 1:
        # Compare exact and thermodynamically-integrated energies:
        plt.figure(2)
        # --- exact energies
        E -= np.interp(0.0, lbda, E)  # difference from bulk
        plt.plot(lbda, E, label="CDFT")
        # --- thermodynamic integration
        dz = z[1] - z[0]
        integrand = (n @ V) * dz
        dlbda = lbda[1] - lbda[0]
        E_TI = hr.trapz(integrand, dlbda)
        E_TI -= np.interp(0.0, lbda, E_TI)  # difference from bulk
        plt.plot(lbda, E_TI, "r+", label="TI")
        plt.axhline(0, color="k", lw=1, ls="dotted")
        plt.axvline(0, color="k", lw=1, ls="dotted")
        plt.legend()
        plt.xlim(lbda.min(), lbda.max())
        plt.xlabel(r"Perturbation strength, $V_0$")
        plt.ylabel(r"Free energy change, $\Delta\Phi$")

    plt.show()


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python plot_data.py <data.h5>")
        exit(1)
    data_file = sys.argv[1]
    plot_data(data_file)


if __name__ == "__main__":
    main()

plt.figure(1)
for i in np.around(np.arange(-10.0,0.1,2),decimals=2):
    with h5py.File(f"water-pl-U{i}-+1.00.h5", "r") as fp:
        r = np.array(fp["r"])
        n = np.array(fp["n"])
        V = np.array(fp["V"])
    plt.plot(r[0::2],n[0::2,1]/(0.033325*2),alpha=1,label=f'$\lambda$ = {i}')
    plt.xlabel("z [$\AA$]",fontsize=11,fontname="Times New Roman")
    plt.ylabel("$n_{H}(z)/n_{bulk}$",fontsize=11,fontname="Times New Roman")
    # plt.legend()
    plt.legend( prop={
            'family': 'Times New Roman', "size": 7, 'stretch': 'normal'})    # plt.title('Classical, -0.4 eV',fontsize=12,fontname="Times New Roman")
    plt.xticks(fontsize=11,fontname="Times New Roman")
    plt.yticks(fontsize=11,fontname="Times New Roman")


plt.xlim([0,12])
plt.ylim([0,4])
figure = plt.gcf()
figure.set_size_inches(3.35, 2.2)
plt.savefig('H_attractive.pdf', dpi=600, bbox_inches='tight')



import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys

if len(sys.argv) < 2:
    print("Usage: python plot.py <file1.h5> [<n_bulk>]")
    exit(1)
    
# add rejection of initial data

filename = sys.argv[1]
n_bulk = float(sys.argv[2]) if (len(sys.argv) > 2) else None

with h5py.File(filename, "r") as fp:
    r = np.array(fp["r"])
    # n = np.array(fp["n"])/n_bulk
    n = np.array(fp["n"])
    V = np.array(fp["V"])

fig, axes = plt.subplots(2, 1, sharex=True, figsize=(6, 8))
axes[0].plot(r, n)
if n_bulk is not None:
    axes[0].axhline(n_bulk, color='k', ls='dotted')
axes[0].set_ylabel("Density")
axes[1].plot(r, V)
axes[1].set_ylabel("Potential")
axes[1].set_xlabel("r")
# plt.show()
plt.savefig(filename[:-4]+".pdf", bbox_inches='tight')
