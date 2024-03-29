"""
  File:   schwarzschild.py
  Author: Thomas Mueller, HdA

    Calculate light rays in Schwarzschild spacetime for an observer approaching
    the black hole.    
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import m4d


rs = 2.0

def calcKsiCrit(r):    
    xi = rs / r;
    if r >= 1.5 * rs:
        ksicrit = np.arcsin(np.sqrt(6.75 * xi * xi * (1.0 - xi)))
    else:
        ksicrit = np.pi - np.arcsin(np.sqrt(6.75 * xi * xi * (1.0 - xi)))
    return ksicrit


obj = m4d.Object()
obj.setMetric("Schwarzschild")

obj.setSolver("GSL_RK4")
obj.setSolverParam("eps_a", 1e-12)
obj.setSolverParam("stepctrl", True)

boxSize = 50.0
boxOffset = 3.0
obj.setSolverParam("lower_bb", -1e12, -boxSize, -boxSize, -boxSize)
obj.setSolverParam("upper_bb", 1e12, boxSize, boxSize, boxSize)

numImages = 300
rMax = 20.0
rMin =  2.1
rStep = (rMax - rMin)/(numImages - 1)


for n in range(numImages):
    fig = plt.figure(figsize=(7,7))
    
    ax = fig.add_subplot(111, title = "Light rays in Schwarzschild spacetime", aspect='equal')    
    plt.axis([-boxSize-boxOffset, boxSize+boxOffset, -boxSize-boxOffset, boxSize+boxOffset])

    ri = rMax - n * rStep
    obj.setInitialPosition(0.0, ri, 0.5 * np.pi, 0.0)
    
    phi = np.linspace(0, 2*np.pi, 200)
    xbg = boxSize * np.cos(phi)
    ybg = boxSize * np.sin(phi)


    for angle in range(0,181,10):
        ksi = np.radians(angle)
        obj.setInitialLocalNullDirection(-1, np.cos(ksi), 0.0, np.sin(ksi))
        obj.calculateGeodesic()
        num = obj.getNumPoints()

        x = []
        y = []

        for i in range(0,num):
            pos = obj.getPosition(i)
            r = pos[1]
            phi = pos[3]
            x.append(r * np.cos(phi))
            y.append(r * np.sin(phi))


        plt.plot(x, y, 'b')
        
        pos = obj.getPosition(num-1)
        if pos[1] > 0.9*boxSize:
            r = boxSize + 2.0
            phi = pos[3]
            ax.text(r * np.cos(phi), r * np.sin(phi), r'${0:.0f}\degree$'.format(angle), fontsize=10, ha='center', va='center')
            
        
    # Milky Way background
    plt.plot(xbg, ybg, 'brown', linewidth=2)
    ax.text(0.55*boxSize, -0.88*boxSize, "Milky Way background", fontsize=10, color='brown')
    
    # black hole
    bh = plt.Circle((0,0), rs, color='black')
    ax.add_patch(bh)
    ax.text(-1.5, -0.7, "BH", fontsize=10, weight='bold', color='white')
        
    # observer with distance annotation
    ax.text(ri, -2, "observer", fontsize=10)
    ax.text(ri, -4, r'$r = {0:.2f}r_s$'.format(ri), fontsize=10)
    ax.text(ri, -6, r'$\xi_c = {0:.2f}\degree$'.format(180.0 - np.degrees(calcKsiCrit(ri))), fontsize=10)
        
    plt.xticks([-40,-20,0,20,40], [r'$-20r_s$', r'$-10r_s$', r'$0r_s$', r'$10r_s$', r'$20r_s$'])
    plt.yticks([-40,-20,0,20,40], [r'$-20r_s$', r'$-10r_s$', r'$0r_s$', r'$10r_s$', r'$20r_s$'])
    
    plt.grid(True)
    
    print("Image {:3d}/{:3d}".format(n, numImages))
    fig.tight_layout(pad=0.3)
    plt.savefig("output/img_{0:03d}.png".format(n), dpi=120, bbox_inches='tight', facecolor='white', pad_inches=0.25)
    plt.close()
