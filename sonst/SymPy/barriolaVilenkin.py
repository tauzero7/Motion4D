"""
   barriolaVilenkin.py 
   
   Start within isympy
   Exec file via  execfile('barriolaVilenkin.py')   
"""

from sympy import *
import m4d

# --- ENTER COORDINATES AS SYMBOLS
t,r,theta,phi = symbols('t r theta phi')

# --- ENTER ADDITIONAL CONSTANTS
k, c = symbols('k c')

# --- SET COORDINATE NAMES AND METRIC COMPONENTS
x_coords = (t,r,theta,phi)
gdd = Matrix((
    (-c**2,0,0,0),
    (0,1,0,0),
    (0,0,k**2*r**2,0),
    (0,0,0,k**2*r**2*sin(theta)**2) ))
    

g_metric     = m4d.Metric(gdd)
christoffels = m4d.Gamma(g_metric,x_coords)

m4d.pprint_christoffels(christoffels)

#m4d.codeprint_metric(g_metric)
#m4d.codeprint_christoffels(christoffels)
#m4d.codeprint_chrisD(christoffels,x_coords)

riem = m4d.Riemann(g_metric,christoffels,x_coords)
m4d.pprint_riemann_down(riem)
#m4d.codeprint_riem(riem)

#ricci = m4d.Ricci(riem)
#m4d.codeprint_ricci(ricci)

#rsc = m4d.RicciScalar(ricci)
#rsc.value()
