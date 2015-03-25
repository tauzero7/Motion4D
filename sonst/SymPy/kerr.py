"""
   kerr.py 
   
   Start within isympy
   Exec file via  execfile('kerr.py')   
"""

from sympy import *
import m4d

# --- ENTER COORDINATES AS SYMBOLS
t,r,theta,phi = symbols('t r theta phi')

# --- ENTER ADDITIONAL CONSTANTS
rs,c,a = symbols('rs c a')

# --- SET COORDINATE NAMES AND METRIC COMPONENTS
x_coords = (t,r,theta,phi)
gdd = Matrix((
    (-(1-rs*r/(r**2+a**2*cos(theta)**2))*c**2,0,0,-rs*a*r*sin(theta)**2/(r**2+a**2*cos(theta)**2)),
    (0, (r**2+a**2*cos(theta)**2)/(r**2-rs*r+a**2),0,0),
    (0,0,r**2+a**2*cos(theta)**2,0),
    (-rs*a*r*sin(theta)**2/(r**2+a**2*cos(theta)**2),0,0,(r**2+a**2+rs*a**2*r*sin(theta)**2/(r**2+a**2*cos(theta)**2))*sin(theta)**2) ))
    

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
