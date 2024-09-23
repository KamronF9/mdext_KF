import os
import h5py
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

# Read H5

filename = '/home/kamron/mdext_KF/examples/water/021_1Dto3Dtest/random1.h5'

# 1D array

with h5py.File(filename, "r") as fp:
    z = np.array(fp["z"])
    # V = np.array(fp["V"])  # only used in plotting
    n = np.array(fp["n"])
    # E = np.array(fp["E"])
    # dE_dn = np.array(fp["dE_dn"])

# print(n.shape) # 16,1,98
# print(n[-1]) 
nTest = n[-1][0]
# print(len(nTest))

# Planar case:
xyGrid = 100
# print(np.broadcast_to(nTest, (xyGrid, len(nTest))))


# Cylindrical case:
# Simple interpolation method
# Make grid
# take the largest z value and set to the length of one side of square
# sin 45 = (a/2) / r_max
# sqrt(2)/2 = a/2 / r_max
# a = r_max*sqrt(2)
a = np.sqrt(2)*z[-1] 
# make grid with 0,0 as the center of the square
nx, ny = (100,100)
x = np.linspace(-a/2, a/2, nx)
y = np.linspace(-a/2, a/2, ny)
gridPts = np.meshgrid(x, y)
xv, yv = gridPts
# plt.plot(xv, yv, marker='o', color='k', linestyle='none')
# plt.savefig('grid.pdf')
# print(np.linalg.norm(gridPts, axis=0).shape) 100,100
r = np.linalg.norm(gridPts, axis=0).flatten() # 1000,

# sys.exit(0)
# r = np.sqrt(x**2+y**2)[0] OK
# print(r.shape) OK
# --- initialize colormap to color by V0
normalize = mpl.colors.Normalize(vmin=r.min(), vmax=r.max())
cmap = mpl.cm.get_cmap("RdBu")
# plt.plot(, , marker='o',c=r, cmap=cmap, linestyle='none')
plt.scatter(xv,yv, marker='o',c=cmap(normalize(r)) )
# --- add colorbar
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=normalize)
sm.set_array([])
plt.colorbar(sm, label=r"radius", ax=plt.gca())
plt.axis('square')
plt.savefig('gridR.pdf')