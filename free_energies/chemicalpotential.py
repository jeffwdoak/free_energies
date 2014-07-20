#!/usr/bin/python

# chemicalpotential.py v0.2 5-22-2012 Jeff Doak jeff.w.doak@gmail.com

# Change log:
# v0.2 - Remade scipt into a class for compatibility with other free energy
#        scripts.
#      - Now using Phase class to input compositions and free energies.

import sys
import numpy as np
from phase import *

def input_energies(file_):
    """
    Read in chemical formulas and total energies (per atom) of a set of
    compounds from the file object file_. Returns a list of phases containing
    composition and energetic data.
    """
    f = open(file_,"r")
    names = []
    stoich = []
    energies = []
    for line in f:
        line = line.split()
        names.append(line[0])
        energies.append(float(line[-1]))
        stoich.append([])
        for i in range(1,len(line)-1):
            stoich[-1].append(float(line[i]))
    phases = []
    for i in range(len(names)):
        phases.append(Phase(stoich[i],energies[i]))
    #return names,stoich,energies
    return phases

class ChemicalPotential():
    """
    Class to calculate the chemical potentials of an n-dimensional system based
    on the equilibrium of n compounds. The compositions of the compounds, as
    well as the free energies of the compounds are expected to be contained in
    instances of the Phase class. Equilibrium can be calculated at any arbitrary
    temperature, assuming that a temperature-dependent free energy is contained
    in the Phase data. ChemicalPotentials class has the following attributes:
        phases - list of compounds in equilibrium. Elements of the phases are
            assumed to be instances of the Phase class
        T - temperature at which to calculate chemical potentials
        n - number of compounds in phases. Must also be the dimensionality of
            composition space
        stoich - nxn numpy array containing the compositions of each compound in
            equilibrium
        energies - nx1 numpy array containing the free energy of each compound
            at the temperature T
        mu - nx1 numpy array containing the chemical potentials of the n
            elements in the compounds
    """
    def __init__(self,phases,T=0):
        self.phases = phases
        self.T = T
        self.n = len(phases)
        self.stoich = np.zeros((self.n,self.n))
        for i in range(len(phases)):
            self.stoich[i] = np.array(phases[i].c)
        energies = []
        for i in range(len(phases)):
            energies.append(phases[i].spline_fit("F")(T))
        self.energies = np.array(energies)
        self.calc_mu()

    def calc_mu(self):
        """
        Calculate the chemical potentials of the n elements in hte n compounds
        in equilibrium.
        """
        self.mu = np.dot(np.linalg.inv(self.stoich),self.energies)
            
if __name__ == "__main__":
    file_ = str(sys.argv[1])
    phases = input_energies(file_)
    chempots = ChemicalPotential(phases)

    print stoich
    print energies
    print chempots.mu
