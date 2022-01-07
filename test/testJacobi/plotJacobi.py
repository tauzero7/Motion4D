"""
  File:   plotJacobi.py
  Author: Thomas Mueller

    Plot geodesic calculated with 'm4dTestJacobi kerr_lens.init > points_lens.dat'

  This file is part of the m4d-library.
"""
import numpy as np
import matplotlib.pyplot as plt

dat = np.genfromtxt("points_lens.dat")

fig = plt.figure()
ax = fig.add_subplot()
ax.set_xlim(0, 1800)
ax.set_ylim(1e-4, 1e3)
ax.set_yscale('log')
ax.plot(dat[:,0], dat[:,8])
ax.grid(True)

plt.savefig("kerr_lens.png", bbox_inches='tight')
#plt.show()
