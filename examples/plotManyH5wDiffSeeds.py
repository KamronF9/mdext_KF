# %%
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import os
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import glob
import matplotlib.ticker as plticker

import matplotlib.gridspec as gridspec
import itertools
# from density_set import DensitySet

def main() -> None:

    particle = 1 # on a 0 index basis
    N_bulk = 0.015


    # Parameter initialization update these if expanding
    # fig_name = 'test'
    fig_name = 'ManySeedsWBMH'  # include traceability here to training version

    # figLabel = ['a', 'b']
    # titles = []  # titles of filename/MD scenario
    # pattern = r'*.h'
    # os.getcwd()

    # Identify all dat files and read them in as data1 and data2
    # First dir sets all the titles up so be sure to have the same dat files in the second dir


    os.chdir(r'/home/kamron/mdext_KF/examples/molten_salt/data9seedtestBMH12NVTH5s/') 

    def GetData(pathAndFilename):

        with h5py.File(pathAndFilename, "r") as fp:
            r = np.array(fp["r"])
            n = np.array(fp["n"])

        return r, n  # list containing 2D layers of each scenario

    

        # Get data

    #     for i, direct in enumerate(directs):

    #         if ax_ind == 0:
    #             #Bottom - repulsive
    #             r,n,errs = GetData(direct, 5.0)  # , potential eV
    #             plt.ylim((0,5))
    #             plt.xlim((0,5.5))
                
    #         else:
    #             #Top - attractive
    #             r,n,errs = GetData(direct, -5.0)
    #             plt.ylim((0,5))
    #             plt.xlim((0,6))
    #             plt.xlabel("z [$\AA$]",fontsize=14)    
    #         density = n[:,particle]/N_bulk
    #         err = errs[:,particle]/N_bulk*2 # set to 2sigma
    #         plt.plot(r, density, label=labels[i])
    #         plt.fill_between(r, density-err, density+err,facecolor='b',alpha=0.5)
    #         plt.text(-0.2, 1.01, f"({figLabel[ax_ind]})", ha="left", va="top",
    #             transform=ax.transAxes, fontsize="large", fontweight="bold")
            
    #         plt.ylabel("$n_{Na}(z)/n_{bulk}$",fontsize=14)

    #         plt.xticks(fontsize=14)
    #         plt.yticks(fontsize=14)

    #     # ax.set_title(titles[j])
    #     ax.grid(True)
    #     # loc = plticker.MultipleLocator(base=2) # this locator puts ticks at regular intervals
    #     # ax.xaxis.set_major_locator(loc)
    #     # loc = plticker.MultipleLocator(base=0.5) # this locator puts ticks at regular intervals
    #     # ax.yaxis.set_major_locator(loc)
    #     # j-=1
    
    # axs[0].legend()



    # fig.savefig(fig_name+'combo.pdf', bbox_inches='tight')
    
    # sys.exit(0)

    # Repulsive

        # pattern = r'*lcurve*'
    pattern = r'*.h'

    plt.figure(figsize=(6,4), dpi=300)
# ----------------
    for pathAndFilename in sorted(glob.iglob(os.path.join(os.getcwd(), pattern))):
        r,n = GetData(pathAndFilename)
        plt.plot(r, n[:,particle]/N_bulk)
        plt.ylim((0,5))
        # plt.legend()
        plt.xlabel("z [$\AA$]")
        plt.ylabel("$n_{Na}(z)/n_{bulk}$")

    plt.savefig(fig_name+'attractWDiffSeeds.pdf', bbox_inches='tight')



    # plt.figure(figsize=(8,6), dpi=300)

    # for i, direct in enumerate(directs):
    #     r,n = GetData(direct, -5.0)
    #     plt.plot(r, n[:,particle]/N_bulk, label=labels[i])
    #     plt.ylim((0,4))
    #     plt.legend()
    #     plt.xlabel("z [$\AA$]",fontsize=11,fontname="Times New Roman")
    #     plt.ylabel("$n_{Na}(z)/n_{bulk}$",fontsize=11,fontname="Times New Roman")

    # plt.savefig(fig_name+'repulse.pdf', bbox_inches='tight', dpi=300)


    # plt.figure(figsize=(8,6), dpi=300)




if __name__ == "__main__":
    main()
