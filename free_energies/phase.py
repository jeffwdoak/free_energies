#!/usr/bin/python

# phase.py v2.4.0 2/16/2012 Jeff Doak jeff.w.doak@gmail.com

# Change log:
# v2.4.0 - Moved to dictionary data structure for storing thermodynamic data
#          by making Phase a dictionary class.
#        - Removed free_energy and _stoopid methods.
#        - Added spline_fit method which generalizes the removed free_energy
#          method, and incorporates the case where _stoopid was needed.
# v2.3.2 - Beginning to work on generalizing spline fitting procedure to any 
#          thermodynamic functions.
#        - Beginning to implement dictionary to store thermodynamic functions

import scipy as sp
from scipy.optimize import leastsq
from scipy.interpolate import UnivariateSpline
import sys

BOLTZCONST = 8.617e-2 #meV/K

class Phase(dict):
    """
        Class containing the free energy of a phase with a given composition.
    Sets of phases can be used to construct a mixing free energy model for a
    region of composition space. Phase class has the following attributes:
        E - zero Kelvin DFT energy (meV/atom)
        c - composition (currently assumed to be pseudobinary between 0 and 1)
        T - temperature (K)
        E_vib - vibrational energy (meV/atom)
        E_ZP - vibrational zero-point energy (meV/atom)
        E_0 - zero Kelvin DFT energy + zero-point energy (meV/atom)
        E_thermal - vibrational energy - zero-point energy (meV/atom)
        S_vib - vibrational entropy (meV/K/atom)
        C_V - heat capacity at constant volume (meV/K/atom)
        F - Helmholtz free energy including vibrations (meV/atom)
    The vibrational attributes and temperature are all arrays with the same
    length. These attributes are stored both as variables within the class, and
    as key - value pairs within the class dictionary.
        In addition to accessing the arrays of thermodynamic data, a cubic
    spline can be fit to any of the data by using the instance method
    spline_fit("key_name")
    where "key_name" is a string containing the name of one of the attributes
    listed above.
        The vibrational thermodynamic data can be read in from an output file of
    the gamma-point frozen phonon code GoBaby written by Vidvus Ozolins at UCLA.
    More file types will be added at some point in the future...
        The vibrational thermodynamic data is also assumed to be for the
    harmonic approximation.
    """
    def __getattr__(self,attr):
        return self[attr]

    def __setattr__(self,attr,val):
        self[attr] = val

    def __init__(self,composition,energy,vibfile=None):
        self.E = energy
        self.c = composition
        # Empty place holders for vibrational thermodynamic data:
        self.T = [0.0]; self.E_vib = [0.0]; self.S_vib = [0.0]; self.C_V = [0.0]
        # Decomposition of vibrational data into zero Kelvin and finite T terms:
        self.E_ZP = 0.0; self.E_thermal = [0.0]; self.E_0 = self.E + self.E_ZP
        self.F = self.E
        # Choose which format vibrational free energies are from, if any:
        if vibfile is not None:
            self.input_GoBaby(open(vibfile,"r"))

    def input_GoBaby(self,file):
        """
        Reads in vibrational free energy data from a GoBaby-formatted file.
        """
        # Assuming for now that file is a file object.
        temp = [];  vibentropy = []; heatcap = []; vibenergy = []
        #file.readline()  # Remove first line containing header information.
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
        # Decompose vibrational energy into zero Kelvin and finite temp terms
        self.E_ZP = self.E_vib[0]
        self.E_thermal = self.E_vib - self.E_ZP
        self.E_0 = self.E + self.E_ZP
        self.F = sp.zeros(len(self.T))
        for i in range(len(self.F)):
            self.F[i] = self.E + self.E_vib[i] - self.T[i]*self.S_vib[i]

    def spline_fit(self,thermo):
        """
        Returns a scipy spline object that is a cubic spline interpolation of
        the thermodynamic data associated with the dictionary key thermo. The
        returned object is a function of one variable, temperature. If the
        thermodynamic data associated with thermo is not temperature dependent,
        the returned function will give the same value back, independent of the
        temperature value given as an input.
        """
        if thermo == "c":
            return self[thermo]
        try:
            len(self[thermo])
            tmp_T = self.T
            tmp_thermo = self[thermo]
        except TypeError:
            tmp_T = sp.linspace(0,2000,201)
            tmp_thermo = sp.ones(201)*self[thermo]
        fit =  UnivariateSpline(tmp_T,tmp_thermo,s=0)
        return fit

def main(argv):
    energyfile = argv[0]
    efile = open(energyfile,"r")
    efile.readline()
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
    #print phase['c']
    #print phase.E
    #print phase['T']
    #print phase.E_vib
    #print phase.F
    #print phase.free_energy
    #print phase.free_energy(phase.T)
    #real = phase.F
    #test = phase.free_energy(phase.T)
    #print real-test
    vib_energy = phase.spline_fit("E_thermal")
    real = vib_energy(phase['T'])
    test = phase["E_thermal"]
    print test

if __name__ == "__main__":
    testPhase(sys.argv[1:])
    #main(sys.argv[1:])
