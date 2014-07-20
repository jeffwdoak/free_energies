#!/usr/bin/python

# binarymixingmodel.py v2.4.1 5/07/2012 Jeff Doak jeff.w.doak@gmail.com

# Change Log:
# v2.4.1 - Made changes to binary_mixing_property to properly deal with
#            multi-dimensional arrays.
# v2.4.0 - Made to work with phase.py v2.4.0.
# v2.3.3 - Added first and second derivative of RK polynomials.
#        - Removed temperature independent RK polynomial function.

import scipy as sp
from scipy.optimize import leastsq
from scipy.interpolate import UnivariateSpline
import sys

BOLTZCONST = 8.617e-2 #meV/K

# Thermodynamic and mixing functions
# Convention: x will correspond to (pseudo-)binary mixing composition ranging
# from 0 to 1, while c will refer to a more general composition metric, possibly
# in a much larger composition space than binary mixing.
def rk_poly(orderX,orderT,vec,x,T):
    """
    Temperature-dependent Redlich-Kister polynomials of any order for 
    two-phase mixing, with a polynomial temperature dependence.
    """
    vec = vec.reshape(orderX+1,orderT+1)
    F_mix = sp.zeros((len(x),len(T)))
    for i in range(len(x)):
        for j in range(len(T)):
            L = sp.zeros(orderX+1)
            for n in range(orderX+1):
                for m in range(orderT+1):
                    L[n] += vec[n,m]*T[j]**m
                F_mix[i,j] += L[n]*(1-2*x[i])**n
            F_mix[i,j] = x[i]*(1-x[i])*F_mix[i,j]
    return F_mix

def d_rk_dx(orderX,orderT,vec,x,T):
    """
    First derivative of the temperature-dependent Redlich Kister polynomials
    with respect to composition.
    """
    vec = vec.reshape(orderX+1,orderT+1)
    dFmix_dx = sp.zeros((len(x),len(T)))
    for i in range(len(x)):
        for j in range(len(T)):
            L = sp.zeros(orderX+1)
            term1 = 0; term2 = 0
            for n in range(orderX+1):
                for m in range(orderT+1):
                    L[n] = vec[n,m]*T[j]**m
                # First term in derivative:
                term1 += L[n]*(1-2*x[i])**(n+1)
            # Second term in derivative:
            for n in range(1,orderX+1):
                term2 += n*L[n]*(1-2*x[i])**(n-1)
            term2 = term2*(-2.)*x[i]*(1-x[i])
            dFmix_dx[i,j] = term1 + term2
    return dFmix_dx

def d2_rk_dx2(orderX,orderT,vec,x,T):
    """
    Second derivative of the temperature-dependent Redlich Kister polynomials
    with respect to composition.
    """
    vec = vec.reshape(orderX+1,orderT+1)
    d2Fmix_dx2 = sp.zeros((len(x),len(T)))
    for i in range(len(x)):
        for j in range(len(T)):
            L = sp.zeros(orderX+1)
            term1 = 0; term2 = 0; term3 = 0
            for n in range(orderX+1):
                for m in range(orderT+1):
                    L[n] = vec[n,m]*T[j]**m
                # First term in derivative:
                term1 += (n+1)*L[n]*(1-2*x[i])**n
            term1 = term1*(-2)
            # Second term in derivative:
            for n in range(1,orderX+1):
                term2 += n*L[n]*(1-2*x[i])**n
            term2 = term2*(-2)
            # Third term in derivative:
            for n in range(2,orderX+1):
                term3 += n*(n-1)*L[n]*(1-2*x[i])**(n-2)
            term3 = term3*4*x[i]*(1-x[i])
            d2Fmix_dx2[i,j] = term1 + term2 + term3
    return d2Fmix_dx2

def binary_mixing_property(property,x,ratio_mixing_atoms=1.0):
    """
    Calculates a thermodnamic property of mixing from an array containing a
    thermodynamic property with units [P]/atom (first index of array
    corresponds to composition) as a function of composition, given by the array
    x. Units of the property of mixing are [P]/mixing_atom, where [P] is the
    unit of the extrinsic property, and the ratio of atoms/mixing_atom is given
    by ratio_mixing_atoms.
    """
    temp = sp.array(property)
    for i in range(1,len(x)-1):
        temp[i] = temp[i]-(1.-x[i])*temp[0]-x[i]*temp[len(x)-1]
    temp[0] = 0.0
    temp[len(x)-1] = 0.0
    temp = temp*ratio_mixing_atoms
    return temp

class BinaryMixingModel:
    """
    Class that creates a mixing free energy model for a pseudo-binary section of
    a composition space.
    """
    def __init__(self,phases,orderX,orderT,num_mixing,key="F"):
        self.orderX = orderX
        self.orderT = orderT
        self.num_mixing = num_mixing
        self.x = []
        self.key = key
        self.thermo = []
        for phase in phases:
            # For now, assume that c = x is a scalar from 0 to 1.
            # This corresponds to mixing in a true binary system.
            # A scheme for determining a pseudo-binary mixing scale from
            # an arbitrary range of composition vectors will need to be 
            # implemented to make this useful.
            self.x.append(phase.c)
            self.thermo.append(phase.spline_fit(self.key))
        self.x = sp.array(self.x)
        if orderT == 0:
            self.T = sp.array([0])
        else:
            self.T = sp.linspace(300,2000,201)  # Let the user pick temp range later!
        # Total free energy of the binary system (without mixing entropy)
        self.F = sp.zeros((len(self.x),len(self.T)))
        for i in range(len(self.x)):
            for j in range(len(self.T)):
                #self.F[i,j] = phases[i].free_energy(self.T[j])
                self.F[i,j] = self.thermo[i](self.T[j])
        # Free energy of mixing in the binary system (without mixing entropy):
        #     Convert units of free energy from meV/atom to meV/mixing_atom. 
        #     This is necessary for adding on an ideal configurational entropy 
        #     term. The composition will also be assumed to be in terms of the 
        #     mixing atoms, not total atoms.
        self.F_mix = binary_mixing_property(self.F,self.x,num_mixing)
        
        #self.init_guess = sp.ones((self.orderX+1)*(self.orderT+1))
        self.init_guess = sp.zeros(((self.orderX+1),(self.orderT+1)))
        self.fit_vec = self.rk_poly_fit(self.init_guess,self.F_mix,self.x,
                           self.T,self.orderX,self.orderT)

    def free_energy(self,x,T):
        """
        Calculate the free energy of mixing in the binary system at an arbitrary
        composition and temperature.
        """
        free = rk_poly(self.orderX,self.orderT,self.fit_vec,x,T)
        free += T*BOLTZCONST*(x*sp.log(x)+(1-x)*sp.log(1-x))
        return free

    def rk_poly_fit(self,init_guess,y,x,T,orderX,orderT):
        """
        Fit binary mixing model to Redlich-Kister polynomials.
        """
        init_guess = init_guess.flatten()
        def residuals(vec,orderX,orderT,y,x,T):
            """
            Residuals of a fit to a binary Redlich-Kister polynomial of a given order.
            """
            err = y - rk_poly(orderX,orderT,vec,x,T)
            return err.flatten()
        fit,success = leastsq(residuals,init_guess,args=(orderX,orderT,y,x,T))
        if not [1,2,3,4].count(success):
            #Error with fitting
            print "There was an error with the solution model fit!"
            sys.exit(1)
        return fit

#This script creates temperature dependent mixing entropy and
#enthalpy coefficients based on a regular or sub-regular solution model.
def main(argv):
    from phase import Phase
    energyfile = argv[0]
    efile = open(energyfile,"r")
    efile.readline()
    energy = float(efile.readline().split()[0])/2.
    inputfile = argv[1]
    phase = Phase(0,energy,inputfile)
    for i in range(len(phase.T)):
        print phase.T[i],2*phase.F[i]/1000.
    sys.exit()

def testBinaryMixingModel(argv):
    """
    Test out the functionality of the new class BinaryMixingModel.
    """
    from phase import Phase
    name_0K = str(argv.pop(0))
    file_0K = open(name_0K,"r")
    file_0K.readline()
    comp = sp.array([])
    energy0 = sp.array([])
    for line in file_0K:
        comp = sp.append(comp,float(line.split()[0]))
        energy0 = sp.append(energy0,float(line.split()[1])/2.)
    phase = []
    for i in range(len(argv)):
        phase.append(Phase(comp[i],energy0[i],str(argv[i])))
    print comp,energy0
    mixing = BinaryMixingModel(phase,1,1,2,"F")
    print mixing.F_mix[1,:]
    print mixing.fit_vec

def testzeroKMixingModel(argv):
    """
    Test out the functionality of the new class BinaryMixingModel.
    """
    from phase import Phase
    name_0K = str(argv.pop(0))
    file_0K = open(name_0K,"r")
    file_0K.readline()
    comp = sp.array([])
    energy0 = sp.array([])
    for line in file_0K:
        comp = sp.append(comp,float(line.split()[0]))
        energy0 = sp.append(energy0,float(line.split()[1])/2.)
    phase = []
    for i in range(len(comp)):
        phase.append(Phase(comp[i],energy0[i]))
    print comp,energy0
    mixing = BinaryMixingModel(phase,1,1,2)
    print mixing.fit_vec


if __name__ == "__main__":
    if len(sys.argv[1:]) == 1:
        testzeroKMixingModel(sys.argv[1:])
    else:
        testBinaryMixingModel(sys.argv[1:])
    #main(sys.argv[1:])
