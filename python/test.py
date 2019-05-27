import m4d

s = m4d.MetricSchwarzschild(1.0)
print(s.getMetricCoeff(0,0))

initPos = m4d.vec4(0.0, 10.0, m4d.M_PI_2, 0.0)

ip = m4d.getTuple(initPos)
print(type(ip))
print(ip)

q = m4d.getThree();
print(q)

#s.calculateMetric(ip)
#s._print()


