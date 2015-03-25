"""
   alcubierreSimple.py 

   Simplified Albubierre Warp metric given by McMonigal et al.
   [PRD85,064024]
   
   Start within isympy
   Exec file via  execfile('alcubierreSimple.py')   
"""

from sympy import *
import m4d

# --- ENTER COORDINATES AS SYMBOLS
t,x,y,z = symbols('t x y z')

# --- ENTER ADDITIONAL CONSTANTS
vs,c = symbols('vs c')

# --- SET COORDINATE NAMES AND METRIC COMPONENTS
x_coords = (t,x,y,z)

U = Function('U')(t,x,y,z)


gdd = Matrix((
    (-c**2+vs**2*U(t,x,y,z)**2,-vs*U(t,x,y,z),0,0),
    (-vs*U(t,x,y,z), 1,0,0),
    (0,0,1,0),
    (0,0,0,1) ))
    

g_metric     = m4d.Metric(gdd)
christoffels = m4d.Gamma(g_metric,x_coords)

m4d.codeprint_metric(g_metric)
m4d.codeprint_christoffels(christoffels)
#m4d.codeprint_chrisD(christoffels,x_coords)

#m4d.pprint_christoffels(christoffels)

#riem = m4d.Riemann(g_metric,christoffels,x_coords)
#m4d.pprint_riemann(riem)
#m4d.codeprint_riem(riem)

#ricci = m4d.Ricci(riem)
#m4d.codeprint_ricci(ricci)

#rsc = m4d.RicciScalar(ricci)
#rsc.value()
