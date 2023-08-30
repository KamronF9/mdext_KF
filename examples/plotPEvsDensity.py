import numpy as np
import pandas as pd
import argparse
import glob
import os
import matplotlib.pyplot as plt

fname = 'densityTotEnPlot'
labels = ('T=0.01','T=0.1','T=1')

# read in csvs, average PE, append to each densities
def get_data():
    data = []
    pattern = r'*csv'
    for pathAndFilename in sorted(glob.iglob(os.path.join(os.getcwd(), pattern))):
        Filename = os.path.basename(pathAndFilename)

        df = pd.read_csv(Filename)
        data.append(np.mean(df.PotEng[100:]))

    data = np.array(data)
    return data

# save mean data
# dfTotEngMean = pd.DataFrame(data, columns=['','TotEngMean'])
# outfile = infile[:-6] + 'csv'
# df.to_csv(outfile)
plt.figure(figsize=(7,5), dpi=300)
x = np.around(np.arange(0.1,1.3,0.1),decimals=1)

for i,folder in enumerate(sorted(glob.glob('*/'))):
    print(folder)
    os.chdir(folder)
    # plot PE vs density
    plt.plot(x, get_data(), rasterized=True, label=labels[i])
    os.chdir('..')
# , label=labels[i]
plt.legend()
plt.ylabel('PE (eV)')
plt.xlabel('n')
# plt.ylim((-5,-2))
# plt.xlim((0,1e4))
# plt.yscale('log')
plt.grid()
# plt.savefig(fname+'.png', bbox_inches='tight')
plt.savefig(fname+'.pdf', bbox_inches='tight')
# plt.show()
# plt.close()