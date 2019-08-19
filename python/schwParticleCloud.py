
import numpy as np
import matplotlib.pyplot as plt
import m4d

obj = m4d.Object()
obj.setMetric("SchwarzschildCart")

obj.setSolver("GSL_RK4")
obj.setSolverParam("eps_a", 1e-8)
obj.setSolverParam("stepctrl", False)

boxSize = 20
obj.setSolverParam("lower_bb", -1e12, -boxSize, -boxSize, -boxSize)
obj.setSolverParam("upper_bb", 1e12, boxSize, boxSize, boxSize)

maxPoints = 5000
obj.setParam("maxNumPoints", maxPoints)


numParticles = 50
data = np.ndarray((numParticles, maxPoints, 4))
for n in range(numParticles):
    print(n)
    r = 0.1 * np.random.rand()
    phi = np.random.rand() * np.pi*2
    obj.setInitialPosition(0.0, 10.0 + r * np.cos(phi), r * np.sin(phi), 0.0)
    obj.setInitialLocalTimeDirection(1, 0.0, 1.0, 0.0, 0.3)
    obj.calculateGeodesic()
    num = obj.getNumPoints()
    
    for k in range(num):
        pos = obj.getPosition(k)
        data[n,k,:] = [pos[0], pos[1], pos[2], pos[3]]


fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
ax.set_xlim([-boxSize, boxSize])
ax.set_ylim([-boxSize, boxSize])

for i in range(0, maxPoints, 100):
    ax.plot(data[:,i,1], data[:,i,2], 'k.', markersize=0.5)

plt.show()
