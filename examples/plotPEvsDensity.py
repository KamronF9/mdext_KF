import numpy as np
import pandas as pd
import argparse
import glob
import os
import matplotlib.pyplot as plt

# read in csvs, average PE, append to each densities
data = []
pattern = r'*csv'
for pathAndFilename in sorted(glob.iglob(os.path.join(os.getcwd(), pattern))):
    Filename = os.path.basename(pathAndFilename)

    df = pd.read_csv(Filename)
    data.append(np.mean(df.PE))
    # plot PE vs density

plt.figure(figsize=(7,5), dpi=300)
    

plt.plot(data[i]['step'], data[i][plotName], rasterized=True, label=labels[i])
plt.legend()
plt.xlabel('Step')
plt.ylabel('RMSE')
# plt.ylim((-5,-2))
plt.xlim((0,1e4))
# plt.xscale('symlog')
plt.yscale('log')
plt.grid()
plt.savefig(fname+'.png', bbox_inches='tight')
plt.savefig(fname+'.pdf', bbox_inches='tight')
plt.show()
plt.close()