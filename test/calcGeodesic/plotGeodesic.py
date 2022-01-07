"""
  File:   plotGeodesic.py
  Author: Thomas Mueller

    Plot geodesic calculated with 'm4dCalcGeodesic kerr.init > points_kerr.dat'

  This file is part of the m4d-library.
"""
import numpy as np
import matplotlib.pyplot as plt

dat = np.genfromtxt("points_kerr.dat")

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_xlim3d(-3, 3)
ax.set_ylim3d(-3, 3)
ax.set_zlim3d(-3, 3)

ax.plot(dat[:,2], dat[:,3], dat[:,4])
plt.savefig("kerr_geod.png", bbox_inches='tight')
