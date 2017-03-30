"""
   schwarzschild.py 
   
   Start within isympy
   Exec file via  execfile('schwarzschild.py')   
"""

from sympy import *
import m4d

# --- ENTER COORDINATES AS SYMBOLS
t,r,theta,phi = symbols('t r theta phi')

# --- ENTER ADDITIONAL CONSTANTS
rs,c = symbols('rs c')

# --- SET COORDINATE NAMES AND METRIC COMPONENTS
x_coords = (t,r,theta,phi)
gdd = Matrix((
    (-(1-rs/r)*c**2,0,0,0),
    (0, 1/(1-rs/r),0,0),
    (0,0,r**2,0),
    (0,0,0,r**2*sin(theta)**2) ))
    

g_metric     = m4d.Metric(gdd)
christoffels = m4d.Gamma(g_metric,x_coords)

m4d.pprint_christoffels(christoffels)
#m4d.codeprint_metric(g_metric)
#m4d.codeprint_christoffels(christoffels)
#m4d.codeprint_chrisD(christoffels,x_coords)

#riem = m4d.Riemann(g_metric,christoffels,x_coords)
#m4d.pprint_riemann(riem)
#m4d.codeprint_riem(riem)

#ricci = m4d.Ricci(riem)
#m4d.codeprint_ricci(ricci)

#rsc = m4d.RicciScalar(ricci)
#rsc.value()
