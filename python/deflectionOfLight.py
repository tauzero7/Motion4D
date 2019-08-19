"""
  File:    deflectionOfLight.py
  Author:  Thomas Mueller
  
  Deflection of light in the Schwarzschild spacetime (cartesian coordinates).
  
  Dependencies:
    - numpy
    - matplotlib
    - m4d library  
"""
import numpy as np
import matplotlib.pyplot as plt
import m4d

obj = m4d.Object()
obj.setMetric("SchwarzschildCart")
obj.setMetricParam("mass", 0.1)

obj.setSolver("GSL_RK4")
obj.setSolverParam("eps_a", 1e-12)
obj.setParam("maxNumPoints", 1200)
obj.setSolverParam("stepctrl", False)

boxSize = 20.0
obj.setSolverParam("lower_bb", -1e12, -boxSize, -boxSize, -boxSize)
obj.setSolverParam("upper_bb", 1e12, boxSize, boxSize, boxSize)


yinit = np.array([2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.534])

def drawGeodesics():
    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(111)
    ax.set_xlim(-12, 12)
    ax.set_ylim(-6, 6)
    ax.set_aspect(1.0)
    ax.grid(True)

    for yi in yinit:        
        # set initial position
        obj.setInitialPosition(0.0, 10.0, yi, 0.0)
        
        # set initial direction (first argument must be integer)
        obj.setInitialLocalNullDirection(-1, -1.0, 0.0, 0.0)
        obj.setSolverParam("stepsize", 0.02)
        
        # calculate geodesic
        obj.calculateGeodesic()

        num = obj.getNumPoints()
        print(yi, num)

        x = []
        y = []
        for i in range(num):
            pos = obj.getPosition(i)
            x.append(pos[1])
            y.append(pos[2])
            
        ax.plot(x, y)
        
    plt.show()


def writeGeodesics():
    for i,yi in enumerate(yinit):
        print(yi)
        # set initial position
        obj.setInitialPosition(0.0, 10.0, yi, 0.0)
        
        # set initial direction (first argument must be integer)
        obj.setInitialLocalNullDirection(-1, -1.0, 0.0, 0.0)
        
        # calculate geodesic
        obj.calculateGeodesic()
        num = obj.getNumPoints()

        data = np.ndarray((num,4))
        for j in range(num):
            pos = obj.getPosition(j)
            data[j,0] = pos[0]
            data[j,1] = pos[1]
            data[j,2] = pos[2]
            data[j,3] = pos[3]
            #print(pos[0], pos[1])
        data.astype(np.float32).tofile("geod_{0}.dat".format(i))
        

drawGeodesics()
#writeGeodesics()        
