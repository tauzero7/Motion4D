"""
  File:   plotSchwGeods.py
  Author: Thomas Mueller

    Plot geodesics calculated with 'm4dTestGeodesic'

  This file is part of the m4d-library.
"""
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot()
ax.set_xlim(-10,10)
ax.set_ylim(-10,10)
ax.set_aspect(1.0)
ax.grid(True)

for i in range(15):
    dat = np.genfromtxt("points_{}.dat".format(i))
    ax.plot(dat[:,2], dat[:,3])

bh = plt.Circle((0,0), 0.5, color='black', zorder=10)
ax.add_artist(bh)

plt.savefig("schw_geods.png", bbox_inches='tight')
#plt.show()
