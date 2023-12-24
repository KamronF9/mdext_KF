import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import glob

# python /home/kamron/mdext_KF/examples/plot1DLJmanyPressures.py

# if len(sys.argv) < 2:
#     print("Usage: python plot.py <file1.h5> [<n_bulk>]")
#     exit(1)
    
plt.figure(figsize=(6, 8))

for i,filename in enumerate(sorted(glob.glob('*h5'))):
# filename = sys.argv[1]
# n_bulk = float(sys.argv[2]) if (len(sys.argv) > 2) else None

    with h5py.File(filename, "r") as fp:
        r = np.array(fp["r"])
        # n = np.array(fp["n"])/n_bulk
        n = np.array(fp["n"])
        # V = np.array(fp["V"])

    label = filename[-6:-3]
    plt.plot(r, n, label='P='+label)
    # if n_bulk is not None:
    #     axes[0].axhline(n_bulk, color='k', ls='dotted')


plt.ylabel("Density")
plt.xlabel("r")
plt.xlim((0,None))
plt.legend()

# plt.show()
plt.savefig('densitiesNT'+".pdf", bbox_inches='tight')
