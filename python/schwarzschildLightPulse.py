"""
  File:   schwarzschildLightPulse.py
  Author: Thomas Mueller, HdA

    Light pulse in Schwarzschild spacetime
"""

import numpy as np
import matplotlib.pyplot as plt
import m4d

obj = m4d.Object()
obj.setMetric("Schwarzschild")

obj.setSolver("GSL_RK4")
obj.setSolverParam("eps_a", 1e-8)
obj.setSolverParam("stepctrl", False)

boxSize = 20.0
obj.setSolverParam("lower_bb", -1e12, -boxSize, -boxSize, -boxSize)
obj.setSolverParam("upper_bb", 1e12, boxSize, boxSize, boxSize)

obj.setInitialPosition(0.0, 10.0, np.pi/2, 0.0)

maxPoints = 1000
obj.setParam("maxNumPoints", maxPoints)

data = np.ndarray((361,maxPoints,2))
for n in range(0,361):
    print(n)
    alpha = 2.0 * np.pi / 360.0 * n;
    obj.setInitialLocalNullDirection(1, -np.cos(alpha), 0.0, np.sin(alpha))
    obj.calculateGeodesic()
    num = obj.getNumPoints()    
    for i in range(num):
        pos = obj.getPosition(i)
        data[n,i] = [pos[1] * np.cos(pos[3]), pos[1] * np.sin(pos[3])]
    

maxLambda = 300
plt.plot(data[:,maxLambda,0],data[:,maxLambda,1],'r.')

maxLambda = 500
plt.plot(data[:,maxLambda,0],data[:,maxLambda,1],'g.')

maxLambda = 700
plt.plot(data[:,maxLambda,0],data[:,maxLambda,1],'b.')

maxLambda = 900
plt.plot(data[:,maxLambda,0],data[:,maxLambda,1],'k*')

plt.show()
#print(data[:,1])
