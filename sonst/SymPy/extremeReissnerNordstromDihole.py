"""
   extremeReissnerNordstromDihole.py 
   
   Start within isympy
   Exec file via  execfile('extremeReissnerNordstromDihole.py')   
"""

from sympy import *
import m4d

# --- ENTER COORDINATES AS SYMBOLS
t,x,y,z = symbols('t x y z')

# --- ENTER ADDITIONAL CONSTANTS
M1,M2 = symbols('M1 M2')

U = Function("U")

# --- SET COORDINATE NAMES AND METRIC COMPONENTS
x_coords = (t,x,y,z)
gdd = Matrix((
    (-1/U(x,y,z)**2,0,0,0),
    (0,U(x,y,z)**2,0,0),
    (0,0,U(x,y,z)**2,0),
    (0,0,0,U(x,y,z)**2) ))
    

g_metric     = m4d.Metric(gdd)
christoffels = m4d.Gamma(g_metric,x_coords)

#m4d.codeprint_metric(g_metric)
#m4d.codeprint_christoffels(christoffels)
#m4d.codeprint_chrisD(christoffels,x_coords)

#m4d.pprint_christoffels(christoffels)

riem = m4d.Riemann(g_metric,christoffels,x_coords)
m4d.pprint_riemann_down(riem)
#m4d.codeprint_riem(riem)

#ricci = m4d.Ricci(riem)
#m4d.codeprint_ricci(ricci)

#rsc = m4d.RicciScalar(ricci)
#rsc.value()
