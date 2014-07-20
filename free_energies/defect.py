#!/usr/bin/python

# defect.py v0.1 5/22/2012 Jeff Doak jeff.w.doak@gmail.com

# Change Log:

import numpy as np
BOLTZCONST = 8.617e-5 #eV/K

class Defect():
    """
        Class containing the energetic, compositional, and structural data about
    a defect structure. Defect class has the following attributes:
        E - raw DFT energy of the defect (eV/u.c.)
        E_ref - raw DFT energy of the perfect (undefected) crystal (eV/u.c.)
        delta_n - nx1 numpy array containing the difference in number of atoms
            between the defected and perfect cell, for each type of atom 
            considered in the system
        delta_e - difference in the number of electrons between the defected and
            perfect cells
        Z - number of ways of introducing the defect into the perfect crystal
        mu - nx1 numpy array containing the atomic chemical potentials of each
            atom type considered in the system (eV)
        mu_e - electron chemical potential of the system (eV)

    Phonons and electrostatic potential corrections may be added later...
    """
    def __init__(self,E,E_ref,delta_n,delta_e,Z=1,mu=None,mu_e=0.0):
        self.E = E
        self.E_ref = E_ref
        self.delta_n = np.array(delta_n)
        self.delta_e = delta_e
        self.Z = Z
        if mu:
            self.mu = np.array(mu)
        else:
            self.mu = np.zeros_like(delta_n)
        self.mu_e = mu_e

    def formation_energy(self):
        """
        Calculate the formation energy of the defect relative to the perfect
        crystal using the atomic and electronic chemical potentials stored in
        the class. Formation energies are returned in eV/defect.
        """
        E_f = self.E - self.E_ref - np.dot(self.delta_n,self.mu)
        E_f -= self.delta_e*(self.mu_e)
        return E_f

    def defect_concentration(self,T):
        """
        Calculate the concentration of defects at the given temperature within
        the dilute limit approximation.
        Note: Entropy associated with different ways of introducing a defect
        into the lattice is currently not accounted for!
        """
        E_f = self.formation_energy()
        c = np.exp(-E_f/(BOLTZCONST*T))
        return c

def input_defect(file_):
    """
    Read in defect information from a file.
    """
    f = open(file_,"r")
    f.readline()  # Discard header line.
    defects = []
    for line in f:
        line = line.split()
        n = len(line) - 5  # Determine the number of atom types to consider.
        E = float(line[1])
        E_ref = float(line[2])
        delta_n = np.zeros(n)
        for i in range(n):
            delta_n[i] = int(line[3+i])
        delta_e = int(line[3+n]); Z = int(line[4+n])
        defects.append(Defect(E,E_ref,delta_n,delta_e,Z))
    f.close()
    return defects

if __name__ == "__main__":
    import sys
    from chemicalpotential import *
    from phase import Phase
    defect_file = str(sys.argv[1])
    chem_pot_file = str(sys.argv[2])
    names,stoich,energies = input_energies(chem_pot_file)
    phases = []
    for i in range(len(names)):
        phases.append(Phase(stoich[i],energies[i]))
    mu = ChemicalPotential(phases).mu
    defects = input_defect(defect_file)
    E_f = []
    for defect in defects:
        defect.mu = mu
        E_f.append(defect.formation_energy())
    print E_f
