"""
   schwarzschildIsotropic.py 
   
   Start within isympy
   Exec file via  execfile('schwarzschildIsotropic.py')   
"""

from sympy import *
import m4d

# --- ENTER COORDINATES AS SYMBOLS
t,x,y,z = symbols('t x y z')

# --- ENTER ADDITIONAL CONSTANTS
rho_s = symbols('rho_s')

# --- SET COORDINATE NAMES AND METRIC COMPONENTS
x_coords = (t,x,y,z)

rho = Function('rho')

gdd = Matrix((
    (-((1-rho_s/rho(x,y,z))/(1+rho_s/rho(x,y,z)))**2,0,0,0),
    (0, (1 + rho_s/rho(x,y,z))**4,0,0),
    (0,0,(1 + rho_s/rho(x,y,z))**4,0),
    (0,0,0,(1 + rho_s/rho(x,y,z))**4)))
    

g_metric     = m4d.Metric(gdd)
christoffels = m4d.Gamma(g_metric,x_coords)

m4d.pprint_christoffels(christoffels)

pprint(factor(simplify(christoffels.ddu(0,0,2)))) 