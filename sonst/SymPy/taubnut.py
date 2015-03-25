"""
   taubnut.py 
   
   Start within isympy
   Exec file via  execfile('taubnut.py')   
"""

from sympy import *
import m4d

# --- ENTER COORDINATES AS SYMBOLS
t,r,theta,phi = symbols('t r theta phi')

# --- ENTER ADDITIONAL CONSTANTS
M,l = symbols('M,l')

# --- SET COORDINATE NAMES AND METRIC COMPONENTS
x_coords = (t,r,theta,phi)

D = Function('D')(r)
S = Function('S')(r)

gdd = Matrix((
    (-D(r)/S(r),0,0,-2*D(r)/S(r)*l*cos(theta)),
    (0, S(r)/D(r),0,0),
    (0,0,S(r),0),
    (-2*D(r)/S(r)*l*cos(theta),0,0,-4*l**2*D(r)/S(r)*cos(theta)**2 + S(r)*sin(theta)**2) ))
    

g_metric     = m4d.Metric(gdd)
christoffels = m4d.Gamma(g_metric,x_coords)


fD = Lambda((r,M,l), r**2-2*M*r-l**2)
fS = Lambda((r,l), r**2+l**2)

#m4d.codeprint_metric(g_metric)

f1 = open('taubnut.txt','w')
for i in [0,1,2,3]:
    for j in [0,1,2,3]:
        for k in [0,1,2,3]:
            #pprint(Eq(Symbol('Chr_%i%i^%i' % (i,j,k)), christoffels.ddu(i,j,k).subs(D(r),fD(r,M,l)).subs(S(r),fS(r,l)).doit().simplify() ))
            print >>f1,"   christoffel[{0}][{1}][{2}] = {3};".format(i,j,k,ccode(christoffels.ddu(i,j,k).subs(D(r),fD(r,M,l)).subs(S(r),fS(r,l)).doit()))



#m4d.pprint_christoffels(christoffels)

#m4d.codeprint_christoffels(christoffels)
#m4d.codeprint_chrisD(christoffels,x_coords)

f1.close()
