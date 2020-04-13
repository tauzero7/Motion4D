"""
   File:   warp.py 
   Author: Thomas Mueller
"""

from sympy import *
import m4d

# --- ENTER COORDINATES AS SYMBOLS
t,x,y,z = symbols('t x y z')

# --- ENTER ADDITIONAL CONSTANTS
#c, vs = symbols('c vs')
c = symbols('c')

vs = Function('vs')(t)
r = Function('r')(t,x,y,z)
f = Function('f')(r)

# --- SET COORDINATE NAMES AND METRIC COMPONENTS
x_coords = (t,x,y,z)
gdd = Matrix((
    (-c**2 + vs**2 * f**2,-vs*f,0,0),
    (-vs*f, 1,0,0),
    (0,0,1,0),
    (0,0,0,1)))
    
def p(i,j,k):
    pprint(Eq(Symbol('Chr_%i%i^%i' % (i,j,k)), simplify(christoffels.ddu(i,j,k))))


g_metric     = m4d.Metric(gdd)
christoffels = m4d.Gamma(g_metric,x_coords)

m4d.pprint_christoffels(christoffels)

