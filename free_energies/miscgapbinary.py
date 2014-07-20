#!/usr/bin/python

# miscgapbinary.py v0.1 1/21/2012 Jeff Doak jeff.w.doak@gmail.com

import scipy as sp
from scipy.optimize import leastsq
import BinaryMixingModel as bmm
#from scipy.interpolate import UnivariateSpline
import sys

BOLTZCONST = 8.617e-2 #meV/K

class MiscGapBinary:
    """
    Class that calculates a pseudo-binary miscibility gap based on an
    analytical solution model of mixing.
    """
    def __init__(self,model):
        self.orderX = model.orderX
        self.orderT = model.orderT
        self.RK_coeff = model.fit_vec

    def system(self,x_vec,T):
        """
        System of equations defining two-phase equilibrium in a pseudo-binary
        alloy.
        """
        f = [bmm.dfmin_dx(self.orderX,self.orderT,x_vec[0],T)
                - bmm.dfmin_dx(self.orderX,self.orderT,x_vec[1],T)]
        f.append((bmm.rk_poly(self.orderX,self.orderT,x_vec[1],T)
                - bmm.rk_poly(self.orderX,self.orderT,x_vec[0],T))
                /(x_vec[1]-x_vec[0])
                - bmm.dfmin_dx(self.orderX,self.orderT,x_vec[0],T))
        return f


