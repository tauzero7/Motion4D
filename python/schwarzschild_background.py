"""
  File:   schwarzschild_background.py
  Author: Thomas Mueller, HdA
  
  Description:
    Calculate light rays in the Schwarzschild spacetime to show the distortion
    of the Milky Way background.
"""

import numpy as np
import matplotlib.pyplot as plt
import m4d

obj = m4d.Object()
obj.setMetric("Schwarzschild")

obj.setSolver("GSL_RK4")
obj.setSolverParam("eps_a", 1e-8)
obj.setSolverParam("stepctrl", True)

background_radius = 20.0
background_offset = 0.5
boxSize = 100.0
observer_pos = 10.0

obj.setSolverParam("lower_bb", -1e12, 0.0, -1e12, -1e12)
obj.setSolverParam("upper_bb", 1e12, boxSize, 1e12, 1e12)

obj.setInitialPosition(0.0, observer_pos, np.pi*0.5, 0.0)

maxPoints = 3000
obj.setParam("maxNumPoints", maxPoints)

initAngles = np.array([
    10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160])


fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111)
ax.set_xlim(-background_radius-4*background_offset, background_radius+4*background_offset)
ax.set_ylim(-background_radius-4*background_offset, background_radius+4*background_offset)
ax.set_aspect(1.0)
#ax.set_xticks(np.arange(0,10,1))

for phi in initAngles:
    print(phi)
    alpha = np.radians(phi)
    obj.setInitialLocalNullDirection(-1, np.cos(alpha), 0.0, np.sin(alpha))
    obj.calculateGeodesic()
    num = obj.getNumPoints()
    x = []
    y = []
    lastAngle = 0.0
    for i in range(num):
        pos = obj.getPosition(i)
        if pos[1] < background_radius+0.5:
            x.append(pos[1] * np.cos(pos[3]))
            y.append(pos[1] * np.sin(pos[3]))
            lastAngle = pos[3]
            
    ax.plot(x,y, 'b', zorder=10)    
    ax.text((background_radius+2*background_offset)*np.cos(lastAngle), 
        (background_radius+2*background_offset)*np.sin(lastAngle),
        "${0:.0f}\degree$".format(phi),
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=14,
        zorder=30)

hide_circle = plt.Circle((0,0), background_radius+0.2, color='#d4aa83', linewidth=6, fill=False, zorder=20)
ax.add_patch(hide_circle)

background_circle = plt.Circle((0,0), background_radius, color='black', linewidth=2, fill=False)
ax.add_patch(background_circle)

bh_disk = plt.Circle((0,0), 2.0, color='black', fill=True)
ax.add_patch(bh_disk)
ax.text(0,0, "BH", color='white',  
    horizontalalignment='center',
    verticalalignment='center',
    fontsize=18)

photon_orbit = plt.Circle((0,0), 3.0, color='black', fill=False, linestyle='dashed')
ax.add_patch(photon_orbit)

observer = plt.Circle((observer_pos,0), 0.2, color='black', fill=True, zorder=20)
ax.add_patch(observer)
ax.text(observer_pos, -0.7, "observer", color='black',
    horizontalalignment='center',
    verticalalignment='center',
    fontsize=14)

plt.grid()
plt.show()

