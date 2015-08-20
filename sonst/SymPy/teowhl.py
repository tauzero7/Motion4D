"""
   teowhl.py 
   
   Start within isympy
   Exec file via  execfile('teowhl.py')   
"""

from sympy import *
import m4d

# --- ENTER COORDINATES AS SYMBOLS
t,l,theta,phi = symbols('t l theta phi')

# --- ENTER ADDITIONAL CONSTANTS
b0,c = symbols('b0 c')

# --- SET COORDINATE NAMES AND METRIC COMPONENTS
x_coords = (t,l,theta,phi)
gdd = Matrix((
    (-c**2. + (l**2. + b0**2.) * sin(theta)**2. * b0**4. * c**2./(4.*(l**2. + b0**2.)**3.), 0, 0, -0.5*(l**2. + b0**2.)*sin(theta)**2. * b0**2. * c/(l**2. + b0**2.)**(3./2.)),
    (0, 1., 0, 0),
    (0, 0, (l**2. + b0**2.), 0),
    (-0.5*(l**2. + b0**2.)*sin(theta)**2. * b0**2. * c/(l**2. + b0**2.)**(3./2.), 0, 0, (l**2. + b0**2.)*sin(theta)**2.)     
    ))
   

g_metric     = m4d.Metric(gdd)
m4d.codeprint_metric(g_metric)

christoffels = m4d.Gamma(g_metric,x_coords)
#m4d.codeprint_christoffels(christoffels)
m4d.pprint_christoffel_ddu(christoffels,0,0,1)