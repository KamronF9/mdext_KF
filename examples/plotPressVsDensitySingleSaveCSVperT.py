import numpy as np
import pandas as pd
import argparse
import glob
import os
import matplotlib.pyplot as plt


# first run bash ../../../parseLammpsLogMdext.sh

fname = 'PressVsVolPlotTest'
# labels = ('T=0.01','T=0.1','T=1')

# read in csvs, average PE, append to each densities
def get_data():
    data = []
    pattern = r'*csv'
    for pathAndFilename in sorted(glob.iglob(os.path.join(os.getcwd(), pattern))):
        Filename = os.path.basename(pathAndFilename)

        df = pd.read_csv(Filename)
        # data.append(np.mean(df.PotEng[100:]))
        data.append(1/np.mean(df.Volume[100:]))
        # data.append(np.mean(df.Press[100:]))

    data = np.array(data)
    return data

# save mean data
# dfTotEngMean = pd.DataFrame(data, columns=['','TotEngMean'])
# outfile = infile[:-6] + 'csv'
# df.to_csv(outfile)
plt.figure(figsize=(7,5), dpi=300)
# x = np.around(np.arange(0.1,1.0+0.1,0.1),decimals=1)
x = np.around(np.arange(0.1,1.0+0.1,0.1),decimals=1) # P

folders = sorted(glob.glob('T*/'))
lenFolders = len(folders)

for i,folder in enumerate(folders):
    label = folder[:-1]
    print(label)
    os.chdir(folder)
    # plot PE vs density
    y = get_data()
    # print(y)
    # Build up list of y save as csv
    # TODO
    if i == 0 or i == lenFolders-1:
        plt.plot(x, y, rasterized=True, label=label)
    else:
        plt.plot(x, y, rasterized=True)
    os.chdir('..')
# , label=labels[i]
plt.legend()
# plt.ylabel('PE')
plt.ylabel('1/V')
plt.xlabel('P')
# plt.ylabel('Pressure') # (p*sigma^3/eps)
# plt.xlabel('Vol')
# plt.ylim((-5,-2))
# plt.xlim((0,1e4))
# plt.yscale('log')
plt.grid()
# plt.savefig(fname+'.png', bbox_inches='tight')
plt.savefig(fname+'.pdf', bbox_inches='tight')
# plt.show()
# plt.close()