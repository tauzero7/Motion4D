
import m4d

obj = m4d.Object()
obj.setMetric("SchwarzschildCart")

obj.setSolver("GSL_RK4")
obj.setSolverParam("eps_a", 1e-8)
obj.setSolverParam("stepctrl", True)

boxSize = 20.0
obj.setSolverParam("lower_bb", -1e12, -boxSize, -boxSize, -boxSize)
obj.setSolverParam("upper_bb", 1e12, boxSize, boxSize, boxSize)

obj.printStatus()

yinit = [ 0.0, 1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, \
7.5, 8.0, 9.0, 10.0 ]

f = open('test.txt', 'w')
for n in range(0, len(yinit)):
    obj.setInitialPosition(0.0, 10.0, yinit[n], 0.0)
    obj.setInitialLocalNullDirection(1, -1.0, 0.0, 0.0)

    obj.calculateGeodesic(2000)
    num = obj.getNumPoints()
    print num
    for i in range(0,num):
        pos = obj.getPosition(i)
        buf = '{0:6.4f} {1:6.4f} {2:6.4f} {3:6.4f}\n'.format(pos[0], pos[1], pos[2], pos[3])
        f.write(buf)
    
    f.write("\n")
        
f.write(buf)        

