import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import os
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as plticker

endRange = 15.0
stepSize = 0.5

# particle = 1 # on a 0 index basis
N_bulk = 1.  # dummy value


# plt.figure(1)
fig, ax = plt.subplots(1,1)
plt.sca(ax)
# --- initialize colormap to color by V0
normalize = mpl.colors.Normalize(vmin=-endRange, vmax=endRange)
# cmap = mpl.cm.get_cmap("RdBu")

colorDict = {
	'red':   ((0.0, 1.0, 1.0), (0.5, 0.6, 0.6), (1.0, 0.0, 0.0)),
	'green': ((0.0, 0.0, 0.0), (0.5, 0.6, 0.6), (1.0, 0.3, 0.3)),
	'blue':  ((0.0, 0.0, 0.0), (0.5, 0.6, 0.6), (1.0, 1.0, 1.0))
}
cmap = LinearSegmentedColormap('RedBlue', colorDict)
# you can then use cmap=cmap kwarg in any plot call

for Ui in np.around(np.arange(0,endRange*2 + stepSize ,stepSize), decimals=1)-endRange:  
    print(f"{Ui:+.1f}")
    
    with h5py.File(f"test-U{Ui:+.1f}.h", "r") as fp:
        r = np.array(fp["r"])
        n = np.array(fp["n"])
        V = np.array(fp["V"])
        if Ui==0.: 
            print(f'N_bulk is {n.mean()=}')
    plt.plot(r,n[:]/(N_bulk),color=cmap(normalize(Ui)), lw=1)
    plt.xlabel("z [$\AA$]",fontsize=14)
    plt.ylabel("$n_{LJ}(z)/n_{bulk}$",fontsize=14)
    # plt.ylabel("$n_{H}(z)/n_{bulk}$",fontsize=11,fontname="Times New Roman")
    
    # plt.legend()
    # plt.legend( prop={
            # 'family': 'Times New Roman', "size": 7, 'stretch': 'normal'})    # plt.title('Classical, -0.4 eV',fontsize=12,fontname="Times New Roman")
    loc = plticker.MultipleLocator(base=2) # this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

# plt.plot(r, V, color="k", lw=1, ls="dashed")  # plot largest shaped V
# plt.xlim([r.min(),r.max()])
plt.xlim([0,5])
plt.ylim([0,1])

# --- add colorbar
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=normalize)
sm.set_array([])
plt.colorbar(sm, label=r"Perturbation strength, $\lambda$")
figure = plt.gcf()
# figure.set_size_inches(3.35, 2.2)
plt.savefig('plotManyExtPots1DLJ.pdf', dpi=600, bbox_inches='tight')

