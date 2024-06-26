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
    fig_name = 'MaxPotComparison_12_12_22ang_AllNVT_nacl_data10_1ANGwAIMD_andnoHpotWseeds'  # include traceability here to training version
    # fig_title = 'Simulation RDFs of NNP Potentials'
    # labels = ['FT', 'Standard NNP','Ext Pot NNP', 'AIMD']  # legend labels
    # directs = [r'/home/kamron/mdext_KF/examples/molten_salt/data6Nima40angBMH/',
    #             r'/home/kamron/mdext_KF/examples/molten_salt/data6Nima40angDPMDreg/',
    #             r'/home/kamron/mdext_KF/examples/molten_salt/data6Nima40angDPMDpert/',
    #             r'/home/kamron/mdext_KF/examples/molten_salt/data7AIMDpert/']
    # labels = [ 'NNP1','NNP0.5' ]  # legend labels
    # labels = ['Classical (Fumi-Tosi)', 'NNP', 'NNP-ext', 'AIMD']  # legend labels
    # labels = ['NNP', 'NNP-ext']  # legend labels
    # labels = ['NNP', 'NNP-ext', 'AIMD']  # legend labels
    labels = ['NNP', 'NNP-extNoHpot', 'AIMD']  # legend labels
    directs = [
        # r'/home/kamron/mdext_KF/examples/molten_salt/data8_12AngBMHNVTH5s1ANG/',
        r'/home/kamron/mdext_KF/examples/molten_salt/data10NVTlongerBoxNNPswSeedsPM/NNP/',
        r'/home/kamron/mdext_KF/examples/molten_salt/data10NVTlongerBoxNNPswSeedsPM/NNPext/',
        # r'/home/kamron/mdext_KF/examples/molten_salt/data10NVTlongerBoxNNPsfromPM/D2ClNaPerturbTrain7r11/',
        r'/home/kamron/mdext_KF/examples/molten_salt/data10NVTlongerBoxAIMDfromPM/',
        
        ]


    figLabel = ['a', 'b']
    # titles = []  # titles of filename/MD scenario
    # pattern = r'*dat'
    # os.getcwd()

    # Identify all dat files and read them in as data1 and data2
    # First dir sets all the titles up so be sure to have the same dat files in the second dir

    # Get titles
    os.chdir(r'/home/kamron/mdext_KF/examples/molten_salt/') 

    def GetData(direct, extPot):

        with h5py.File(direct+f"Alltest-U{extPot:+.1f}.h", "r") as fp:
            r = np.array(fp["r"])
            n = np.array(fp["n"])
            try:
                errs = np.array(fp["std"])
            except:
                errs = np.zeros_like(n)
            # V = np.array(fp["V"])

        return r, n, errs  # list containing 2D layers of each scenario

    
    fig, axs = plt.subplots(2, 1, figsize=(5, 7), dpi=300, sharex=True)
    plt.subplots_adjust(hspace=0.13)

    for ax_ind, ax in enumerate(axs.flat):  # order somehow manages to be lq hp crystal so reorder it
        #    print(ax)
        #     ax.subplot(2,3,i+1, sharex=ax1, sharey=ax1)
        # Get data

        for i, direct in enumerate(directs):
            plt.sca(ax)
            if ax_ind == 0:
                #Bottom - repulsive
                r,n,errs = GetData(direct, 5.0)  # , potential eV
                plt.ylim((0,5))
                # plt.xlim((0,5.5))
                
            else:
                #Top - attractive
                r,n,errs = GetData(direct, -5.0)
                plt.ylim((0,5))
                # plt.xlim((0,6))
                plt.xlabel("z [$\AA$]",fontsize=14)    
            density = n[:,particle]/N_bulk
            err = errs[:,particle]/N_bulk*10 # set to 2sigma XXXXXXXXX HACK
            plt.plot(r, density, label=labels[i])
            
            #remove error for now from AIMD
            if i == 0:
                plt.fill_between(r, density-err, density+err,facecolor='blue',alpha=0.2)
            if i == 1:
                plt.fill_between(r, density-err, density+err,facecolor='orange',alpha=0.2)

            plt.text(-0.2, 1.01, f"({figLabel[ax_ind]})", ha="left", va="top",
                transform=ax.transAxes, fontsize="large", fontweight="bold")
            
            plt.ylabel("$n_{Na}(z)/n_{bulk}$",fontsize=14)

            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)

        # ax.set_title(titles[j])
        ax.grid(True)
        # loc = plticker.MultipleLocator(base=2) # this locator puts ticks at regular intervals
        # ax.xaxis.set_major_locator(loc)
        # loc = plticker.MultipleLocator(base=0.5) # this locator puts ticks at regular intervals
        # ax.yaxis.set_major_locator(loc)
        # j-=1
    
    axs[0].legend()



    fig.savefig(fig_name+'combo.pdf', bbox_inches='tight')
    
    # sys.exit(0)

    # Repulsive
    # plt.figure(figsize=(8,6), dpi=300)
#----------------
    # for i, direct in enumerate(directs):
    #     r,n = GetData(direct, 5.0)
    #     plt.plot(r, n[:,particle]/N_bulk, label=labels[i])
    #     plt.ylim((0,4))
    #     plt.legend()
    #     plt.xlabel("z [$\AA$]",fontsize=11,fontname="Times New Roman")
    #     plt.ylabel("$n_{Na}(z)/n_{bulk}$",fontsize=11,fontname="Times New Roman")

    # plt.savefig(fig_name+'attract.pdf', bbox_inches='tight')



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
