# module load venv/mdextKF

import os
import h5py
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys


# # Create a meshgrid
# x = np.linspace(0, 10, 3)  #was 0 10 2
# y = np.linspace(0, 1, 3)
# X, Y = np.meshgrid(x, y)

# # print(X)
# # print(Y)
# # [[ 0. 10.]
# #  [ 0. 10.]
# #  [ 0. 10.]]
# # [[0.  0. ]
# #  [0.5 0.5]
# #  [1.  1. ]]

# # Access a specific row (e.g., the 1st row)
# row = X[0, :]
# print("Row from X:", row)

# # Access a specific column (e.g., the 2nd column)
# column = Y[:, 1]
# print("Column from Y:", column)

# plt.plot(X[0,:], Y[:,0]) 
# plt.savefig('_xySimpleline.pdf')

# sys.exit(1)


# Read H5

filename = '/home/kamron/mdext_KF/examples/water/021_1Dto3Dtest/random1.h5'

# 1D array

with h5py.File(filename, "r") as fp:
    z = np.array(fp["z"])
    V = np.array(fp["V"])  # only used in plotting
    n = np.array(fp["n"])
    E = np.array(fp["E"])
    dE_dn = np.array(fp["dE_dn"])

# print(n.shape) # 16,1,98
# print(n[-1]) 
nTest = n[-1][0] # largest lambda value to test
# print(len(nTest))

# Planar case:
xyPts = 100
# print(np.broadcast_to(nTest, (xyPts, len(nTest))))
# set scale same as z
a = z[-1]
xyGrid = np.linspace(-a/2, a/2, xyPts)
# testY = np.linspace(-a, a, xyPts)
# print(xyGrid)
# plt.plot(z,nTest)
# plt.savefig('a.pdf')

# gridPts = np.meshgrid(xyGrid, testY, z)
gridPts = np.meshgrid(xyGrid, xyGrid, z)
xv, yv, zv = gridPts
gridStack = np.stack(gridPts)
print(gridStack.shape)
# print(xv.shape) #100,100,98
nv = np.broadcast_to(nTest, (xv.shape))

# i=10
# print(xv[i,i,i],yv[i,i,i],zv[i,i,i],nv[i,i,i])

# plt.plot(xv[0,:,0], yv[:,0,0]) # xy line for 0 row/column entry and 0 z
# plt.savefig('_xyline.pdf')

# plt.plot(zv[0,0,:], nv[0,0,:]) 
# print('xy values', xv[10,10,0],yv[99,10,0]) # -3.87, 4.85 OK
# print('z value', zv[10,10,20]) # 2 OK
# print(z) # 0-9.7 OK
# plt.plot(zv[10,10,:], nv[10,10,:]) 
# plt.savefig('_z_nvals.pdf') OK shows the z profile of n
# plt.plot(xv[10,:,10], nv[10,:,10]) 
# plt.savefig('_x_nvals.pdf') OK flat line
# plt.scatter(xv[0,:,:], zv[0,:,:]) # x - -4.85 to 4.85 and y 0 to 9.7 OK
# plt.scatter(yv[:,0,:], zv[0,:,:]) # x - -9.7 to 9.7 with testY grid and y 0 to 9.7 OK
plt.contourf(xv[0,:,:], zv[0,:,:], nv[0,:,:]) 
plt.axis('scaled')
plt.colorbar()
plt.show()
plt.savefig('_xz_nvals.pdf')



'''

out_file = '/home/kamron/mdext_KF/examples/water/021_1Dto3Dtest/test3D.h5'
# Write hdf5 file (back in original units to keep things O(1)):
with h5py.File(out_file, "w") as fp:
    fp["grid"] = gridStack
    fp["V"] = V  # only used in plotting
    fp["n"] = 
    fp["E"] = E_ex / V_unit
    fp["dE_dn"] = V_ex[:, None] / V_unit
# print(
#     f"Wrote {out_file} with grid length {L/r_unit:.2f} A"
#     f" and spacing {dr/r_unit:.2f} A."
# )

sys.exit(0)


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


# random.h5:
# HDF5 "/home/kamron/mdext_KF/examples/water/021_1Dto3Dtest/random1.h5" {
# GROUP "/" {
#    DATASET "E" {
#       DATATYPE  H5T_IEEE_F64LE
#       DATASPACE  SIMPLE { ( 16 ) / ( 16 ) }
#    }
#    DATASET "V" {
#       DATATYPE  H5T_IEEE_F64LE
#       DATASPACE  SIMPLE { ( 1, 98 ) / ( 1, 98 ) }
#    }
#    DATASET "dE_dn" {
#       DATATYPE  H5T_IEEE_F64LE
#       DATASPACE  SIMPLE { ( 16, 1, 98 ) / ( 16, 1, 98 ) }
#    }
#    DATASET "n" {
#       DATATYPE  H5T_IEEE_F64LE
#       DATASPACE  SIMPLE { ( 16, 1, 98 ) / ( 16, 1, 98 ) }
#    }
#    DATASET "z" {
#       DATATYPE  H5T_IEEE_F64LE
#       DATASPACE  SIMPLE { ( 98 ) / ( 98 ) }
#    }
# }
# }

'''