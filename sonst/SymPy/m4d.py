"""
   script is partially based on the code by 
   https://code.google.com/p/pythonxy/source/browse/src/python/sympy/DOC/examples/advanced/relativity.py?spec=svn.xy-27.b2b5aae1e3df18cf9b4025c528b658faca5b813f&repo=xy-27&r=b2b5aae1e3df18cf9b4025c528b658faca5b813f
  
   Reload module via:  reload(m4d)
"""

import sys
from sympy import *


class Metric(object):
    """ 
        Turn matrix into upper and lower metric
    """
    def __init__(self,m):
        self.gdd = m
        self.guu = m.inv()
        
    def __str__(self):
        return "g_dd = \n" + str(self.gdd)
        
    def dd(self,i,j):
        return self.gdd[i,j]
        
    def uu(self,i,j):
        return self.guu[i,j]
        

class Gamma(object):
    """
       Calculate Christoffel Gamma_ij^k symbols of metric g
    """
    def __init__(self,g,x):
        self.g = g
        self.x = x
        
    def ddu(self,i,j,k):
        g = self.g
        x = self.x
        chr = 0
        for m in [0,1,2,3]:
            chr += g.uu(k,m)/2 * (g.dd(m,i).diff(x[j]) \
                   + g.dd(m,j).diff(x[i]) - g.dd(i,j).diff(x[m]))
        #return chr.simplify()  
        return chr


class Riemann(object):
    """
        Calculate Riemann tensor R^mu_nu,rho,sigma        
    """
    def __init__(self,g,G,x):
        self.g = g
        self.G = G
        self.x = x
        
    def uddd(self,mu,nu,rho,sigma):
        G = self.G
        x = self.x
        R = G.ddu(nu,sigma,mu).diff(x[rho]) - G.ddu(nu,rho,mu).diff(x[sigma])
        for lam in [0,1,2,3]:
            R += G.ddu(rho,lam,mu)*G.ddu(nu,sigma,lam) - \
                 G.ddu(sigma,lam,mu)*G.ddu(nu,rho,lam)
        return R.simplify()
        
    def dddd(self,mu,nu,rho,sigma):
        g = self.g
        R = 0
        for lam in [0,1,2,3]:
            R += g.dd(mu,lam)*self.uddd(lam,nu,rho,sigma)
        return R.simplify()
        

class Ricci(object):
    """
        Calculate Ricci tensor from Riemann tensor
    """
    def __init__(self,R):
        self.R = R
        
    def dd(self,mu,nu):
        R = self.R
        Ric = 0
        for rho in [0,1,2,3]:
            Ric += R.uddd(rho,mu,rho,nu)
        return Ric.simplify()
        
        
class RicciScalar(object):
    """
        Calculate Ricci scalar from Ricci tensor
    """
    def __init__(self,Ric,g):
        self.Ric = Ric
        self.g   = g
        
    def value(self):
        Ric = self.Ric
        g   = self.g
        RS = 0
        for mu in [0,1,2,3]:
            for nu in [0,1,2,3]:
                RS += g.uu(mu,nu)*Ric.dd(mu,nu)
        return RS.simplify()
        

def pprint_christoffel_ddu(Gamma,i,j,k):
    pprint(Eq(Symbol('Chr_%i%i^%i' % (i,j,k)), Gamma.ddu(i,j,k)))


def pprint_christoffels(Gamma):
    for i in [0,1,2,3]:
        for j in [0,1,2,3]:
            for k in [0,1,2,3]:
                if (Gamma.ddu(i,j,k)!=0):
                    pprint_christoffel_ddu(Gamma,i,j,k)

def pprint_riemann(Riemann):
    for i in [0,1,2,3]:
        for j in [0,1,2,3]:
            for k in [0,1,2,3]:
                for m in [0,1,2,3]:
                    if (Riemann.uddd(i,j,k,m)!=0):
                        pprint(Eq(Symbol('R^%i_%i%i%i' % (i,j,k,m)), Riemann.uddd(i,j,k,m)))

def pprint_riemann_down(Riemann):
    for i in [0,1,2,3]:
        for j in [0,1,2,3]:
            for k in [0,1,2,3]:
                for m in [0,1,2,3]:
                    if (Riemann.uddd(i,j,k,m)!=0):
                        pprint(Eq(Symbol('R_%i%i%i%i' % (i,j,k,m)), Riemann.dddd(i,j,k,m)))

def pprint_ricci(Ricci):
    for i in [0,1,2,3]:
        for j in [0,1,2,3]:
            if (Ricci.dd(i,j)!=0):
                pprint(Eq(Symbol('R_%i%i' % (i,j)), Ricci.dd(i,j)))

def codeprint_metric(g,f=sys.stdout):
    for i in [0,1,2,3]:
        for j in [0,1,2,3]:
            #print >>f, "g_compts[{0}][{1}] = {2};".format(i,j,ccode(g.dd(i,j)))
            print("g_compts[{0}][{1}] = {2};".format(i,j,ccode(g.dd(i,j))))


def codeprint_christoffels(Gamma,f=sys.stdout):
    for i in [0,1,2,3]:
        for j in [0,1,2,3]:
            for k in [0,1,2,3]:
                #print >>f, "christoffel[{0}][{1}][{2}] = {3};".format(i,j,k,ccode(Gamma.ddu(i,j,k)))
                print("christoffel[{0}][{1}][{2}] = {3};".format(i,j,k,ccode(Gamma.ddu(i,j,k))))

        
def codeprint_chrisD(Gamma,X,f=sys.stdout):
    for i in [0,1,2,3]:
        for j in [0,1,2,3]:
            for k in [0,1,2,3]:
                for m in [0,1,2,3]:
                    print >>f, "chrisD[{0}][{1}][{2}][{3}] = {4};".format(i,j,k,m,ccode(Gamma.ddu(i,j,k).diff(X[m]).simplify()))


def codeprint_riem(Riemann,f=sys.stdout):
    for i in [0,1,2,3]:
        for j in [0,1,2,3]:
            for k in [0,1,2,3]:
                for m in [0,1,2,3]:
                    print >>f, "riem[{0}][{1}][{2}][{3}] = {4};".format(i,j,k,m,ccode(Riemann.uddd(i,j,k,m)))


def codeprint_ricci(Ricci,f=sys.stdout):
    for i in [0,1,2,3]:
        for j in [0,1,2,3]:
            print >>f, "ricci[{0}][{1}] = {2}".format(i,j,ccode(Ricci.dd(i,j)))
            
