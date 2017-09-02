"""
  Light pulse in Schwarzschild spacetime
"""

import numpy as np
import matplotlib.pyplot as plt
import m4d

obj = m4d.Object()
obj.setMetric("SchwarzschildCart")

obj.setSolver("GSL_RK4")
obj.setSolverParam("eps_a", 1e-8)
obj.setSolverParam("stepctrl", False)

boxSize = 20.0
obj.setSolverParam("lower_bb", -1e12, -boxSize, -boxSize, -boxSize)
obj.setSolverParam("upper_bb", 1e12, boxSize, boxSize, boxSize)

obj.setInitialPosition(0.0, 10.0, 0.0, 0.0)

maxPoints = 1000
obj.setParam("maxNumPoints", maxPoints)

data = np.ndarray((361,maxPoints,4))
for n in range(0,361):
    print(i)
    alpha = 2.0 * np.pi / 360.0 * n;
    obj.setInitialLocalNullDirection(1, -np.cos(alpha), np.sin(alpha), 0.0)
    obj.calculateGeodesic()
    num = obj.getNumPoints()    
    for i in range(num):
        pos = obj.getPosition(i)
        data[n,i] = [pos[0], pos[1], pos[2], pos[3]]
    

maxLambda = 300
plt.plot(data[:,maxLambda,1],data[:,maxLambda,2],'r.')

maxLambda = 500
plt.plot(data[:,maxLambda,1],data[:,maxLambda,2],'g.')

maxLambda = 700
plt.plot(data[:,maxLambda,1],data[:,maxLambda,2],'b.')

maxLambda = 900
plt.plot(data[:,maxLambda,1],data[:,maxLambda,2],'k*')

plt.show()
#print(data[:,1])
