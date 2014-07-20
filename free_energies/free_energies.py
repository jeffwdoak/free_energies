#!/usr/bin/python

# free_energies.py v2.2 12/21/2011 Jeff Doak jeff.w.doak@gmail.com

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
    if orderT == 0:
        F_mix = rk_poly0(orderX,vec,x)
        return F_mix
    #L = sp.zeros(orderX+1)
    #L = sp.zeros((orderX+1,len(T)))
    #F_mix = 0
    #F_mix = sp.zeros(len(T))
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
    #for i in range(orderX+1):
    #    for j in range(orderT+1):
    #        L[i] += vec[i,j]*T**j
    #    print L
    #    print (1-2*x)**i
    #    print (1-2*x)**i*L[i]
    #    F_mix += L[i]*(1-2*x)**i
    #F_mix = x*(1-x)*F_mix
    #print F_mix
    #print len(F_mix),len(F_mix[0])
    return F_mix

def rk_poly0(order,vec,x):
    """
    Temperature-indepenent Redlich-Kister polynomials of any order for
    two-phase mixing.
    """
    L_tot = 0
    for i in range(order+1):
        L_tot += vec[i]*(1-2*x)**i
    F_mix = x*(1-x)*L_tot
    return F_mix

def binary_mixing_property(property,x,ratio_mixing_atoms):
    """
    Calculates a thermodnamic property of mixing from an array containing a
    thermodynamic property with units [P]/atom (first index of array
    corresponds to composition) as a function of composition, given by the array
    x. Units of the property of mixing are [P]/mixing_atom, where [P] is the
    unit of the extrinsic property, and the ratio of atoms/mixing_atom is given
    by ratio_mixing_atoms.
    """
    for i in range(1,len(x)-1):
        property[i] = property[i]-(1.-x[i])*property[0]-x[i]*property[len(x)-1]
    property[0] = 0.0
    property[len(x)-1] = 0.0
    property = property*ratio_mixing_atoms
    return property

class Phase:
    """
    Class containing the free energy of a phase with a given composition. Sets
    of phases can be used to construct a mixing free energy model for a region
    of composition space. Phase class has the following atributes:
        E - zero kelvin energy (meV/atom)
        c - composition (currently assumed to be pseudobinary between 0 and 1)
        T - temperature (K)
        E_vib - vibrational energy (meV/atom)
        S_vib - vibrational entropy (meV/K/atom)
        C_V - heat capacity at constant volume (meV/K/atom)
        F - Helmholtz free energy including vibrations (meV/atom)
    The vibrational atributes and temperature are all arrays with the same
    length. In addition to atributes, there are some instance methods that can
    be called by the user. These include:
        free_energy(T) - cubic spline free energy model for the phase which
            returns a free energy (meV/atom) for any input temperature (K).
    """
    def __init__(self,composition,energy,vibfile=None):
        self.E = energy
        self.c = composition
        #self.T = None; self.E_vib = None; self.S_vib = None; self.C_V = None
        self.T = [0.0]; self.E_vib = [0.0]; self.S_vib = [0.0]; self.C_V = [0.0]
        if vibfile is not None:
            self.input_GoBaby(open(vibfile,"r"))
            self.F()
            self.free_energy = UnivariateSpline(self.T,self.F,s=0)
        else:
            self.F()
            self.free_energy = self._stoopid

    def input_GoBaby(self,file):
        """
        Reads in vibrational free energy data from a GoBaby-formatted file.
        """
        # Assuming for now that file is a file object.
        temp = [];  vibentropy = []; heatcap = []; vibenergy = []
        file.readline()  # Remove first line containing header information.
        for line in file:
            line = line.split()
            if float(line[0]) == 1.0:
                temp.append(0.0)
                vibentropy.append(float(line[1]))
                heatcap.append(float(line[2]))
                vibenergy.append(float(line[3]))
            else:
                temp.append(float(line[0]))
                vibentropy.append(float(line[1]))
                heatcap.append(float(line[2]))
                vibenergy.append(float(line[3]))
        file.close()
        self.T = sp.array(temp)
        self.E_vib = sp.array(vibenergy)
        self.S_vib = sp.array(vibentropy)*BOLTZCONST
        self.C_V = sp.array(heatcap)*BOLTZCONST

    def _stoopid(self,temp):
        """
        Kludgy method for creating a "function" that returns the DFT total
        energy given any temperature. This is to be used in place of a cubic
        spline interpolation of vibrational data when none exists.
        """
        return self.F[0]

    def F(self):
        self.F = sp.zeros(len(self.T))
        for i in range(len(self.F)):
            self.F[i] = self.E + self.E_vib[i] - self.T[i]*self.S_vib[i]

class BinaryMixingModel:
    """
    Class that creates a mixing free energy model for a pseudo-binary section of
    a composition space.
    """
    def __init__(self,phases,orderX,orderT,num_mixing):
        self.orderX = orderX
        self.orderT = orderT
        self.num_mixing = num_mixing
        self.x = []
        for phase in phases:
            self.x.append(phase.c)  # For now, assume c is scalar from 0 to 1.
        self.x = sp.array(self.x)
        self.T = sp.linspace(0,2000,201)  # Let the user pick temp range later!
        # Total free energy of the binary system (without mixing entropy)
        #self.F = sp.zeros((len(self.x),len(self.T)))
        #for i in range(len(self.x)):
        #    for j in range(len(self.T)):
        #        self.F[i,j] = phase[i].free_energy(self.T[j])
        # Free energy of mixing in the binary system (without mixing entropy)
        self.F_mix = sp.zeros((len(self.x),len(self.T)))
        for j in range(len(self.T)):
            for i in range(1,len(self.x)-1):
                self.F_mix[i,j] = (phases[i].free_energy(self.T[j]) -
                        (1.-self.x[i])*phases[0].free_energy(self.T[j]) -
                        self.x[i]*phases[len(self.x)-1].free_energy(self.T[j]))
            self.F_mix[0,j] = 0.0
            self.F_mix[len(self.x)-1,j] = 0.0
        # Convert units of free energy from meV/atom to meV/mixing_atom. This is
        # necessary for adding on an ideal configurational entropy term. The
        # composition will also be assumed to be in terms of the mixing atoms,
        # not total atoms.
        self.F_mix = self.F_mix*self.num_mixing
        
        #self.init_guess = sp.ones((self.orderX+1)*(self.orderT+1))
        self.init_guess = sp.ones(((self.orderX+1),(self.orderT+1)))
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
            #print sp.shape(err.flatten())
            print err.flatten()
            return err.flatten()
        output = leastsq(residuals,init_guess,args=(orderX,orderT,y,x,T))
        return output

#This script creates temperature dependent mixing entropy and
#enthalpy coefficients based on a regular or sub-regular solution model.
def parseargs(argv):
    """
    Function to parse input arguments to script.
    """
    mixingatom = 1
    files = []
    while len(argv) > 0:
        temp = argv.pop(0)
        if str(temp) == "-m" or str(temp) == "--mixing":
            mixingatom = int(argv.pop(0))
        else:
            files.append(temp)
    return mixingatom,files

def main(argv):
    energyfile = argv[0]
    efile = open(energyfile,"r")
    energy = float(efile.readline().split()[0])/2.
    inputfile = argv[1]
    phase = Phase(0,energy,inputfile)
    for i in range(len(phase.T)):
        print phase.T[i],2*phase.F[i]/1000.
    sys.exit()

def testPhase(argv):
    """
    Test out the functionality and validity of the new class Phase.
    """
    if len(argv) > 0:
        inputfile = argv[0]
        phase = Phase(0,0,inputfile)
    else:
        phase = Phase(0,3)
    print phase.c
    print phase.E
    print phase.T
    print phase.E_vib
    print phase.F
    print phase.free_energy
    print phase.free_energy(phase.T)
    real = phase.F
    test = phase.free_energy(phase.T)
    print real-test

def testBinaryMixingModel(argv):
    """
    Test out the functionality of the new class BinaryMixingModel.
    """
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
    mixing = BinaryMixingModel(phase,1,0,2)
    print mixing


if __name__ == "__main__":
    #testBinaryMixingModel(sys.argv[1:])
    #testPhase(sys.argv[1:])
    main(sys.argv[1:])
