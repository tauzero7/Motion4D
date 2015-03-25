"""
   rotDihole.py 
   
   Start within isympy
   Exec file via  execfile('rotDihole.py')   
"""

from sympy import *
import m4d

# --- ENTER COORDINATES AS SYMBOLS
t,x,y,z = symbols('t x y z')

# --- ENTER ADDITIONAL CONSTANTS
M1,M2,omega = symbols('M1 M2 omega')

# --- SET COORDINATE NAMES AND METRIC COMPONENTS
x_coords = (t,x,y,z)

#fU = Lambda((M1,M2,t,x,y,z), 1+M1/sqrt(x**2+y**2+(z-1)**2) + M2/sqrt(x**2+y**2+(z+1)**2));
fU = Lambda((M1,M2,omega,t,x,y,z), \
        1 + M1/sqrt(x**2+(y-sin(omega*t))**2+(z-cos(omega*t))**2) \
          + M2/sqrt(x**2+(y+sin(omega*t))**2+(z+cos(omega*t))**2));

U = Function('U')(t,x,y,z)

gdd = Matrix((
    (-1/U(t,x,y,z)**2,0,0,0),
    (0,U(t,x,y,z)**2,0,0),
    (0,0,U(t,x,y,z)**2,0),
    (0,0,0,U(t,x,y,z)**2) ))
    

g_metric     = m4d.Metric(gdd)
christoffels = m4d.Gamma(g_metric,x_coords)

#m4d.pprint_christoffels(christoffels)


f1 = open('dihole.txt','w')
for i in [0,1,2,3]:
    for j in [0,1,2,3]:
        print >>f1, "    g_compts[{0}][{1}] = {2};".format(i,j,ccode(g_metric.dd(i,j)))
print >>f1,"\n"

for i in [0,1,2,3]:
    for j in [0,1,2,3]:
        for k in [0,1,2,3]:
            #print >>f1, "    christoffel[{0}][{1}][{2}] = {3};".format(i,j,k,ccode(christoffels.ddu(i,j,k).subs(U(t,x,y,z),fU(M1,M2,t,x,y,z)).doit()))
            print >>f1, "    christoffel[{0}][{1}][{2}] = {3};".format(i,j,k,ccode(christoffels.ddu(i,j,k)))

print >>f1, "\n"
print >>f1, "  U = {0};".format(ccode(fU(M1,M2,omega,t,x,y,z)))
print >>f1, "\n"
print >>f1, "  dUdt = {0};".format(ccode(fU(M1,M2,omega,t,x,y,z).diff(t)))
print >>f1, "  dUdx = {0};".format(ccode(fU(M1,M2,omega,t,x,y,z).diff(x)))
print >>f1, "  dUdy = {0};".format(ccode(fU(M1,M2,omega,t,x,y,z).diff(y)))
print >>f1, "  dUdz = {0};".format(ccode(fU(M1,M2,omega,t,x,y,z).diff(z)))
f1.close()

