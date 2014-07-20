#!/usr/bin/python

# charged_defectspy v0.1 5/22/2012 Jeff Doak jeff.w.doak@gmail.com

# Change Log:

# Python script to calculate the formation energy and equilibrium concentrations
# of charged defects, conduction electrons, and valence holes, as a function of
# atomic chemical potential and temperature.
# Script makes use of the ElectronicDOS and Defect classes to deal with
# electrons and atomic defects, respectively. The ChemicalPotential class is
# also used to find the chemical potentials to use in calculating defect
# formation energies.

# Program flow:
# - Read in data
#   - electronic dos
#   - file with atomic defects
#   - file with compound energetics for chemical potentials
# - Calculate atomic chemical potentials
# - Calculate electonic chemical potentials

import numpy as np
import sys
from defect import *
from atat_electron import *  # Change this when I move atat_electron to electronicdos
from chemicalpotential import *
from scipy.optimize import brentq
BOLTZCONST = 8.617e-5 #eV/K

def charge_neutrality(mu_e,T,defects,doscar):
    """
    Calculates the condition of charge neutrality between charged defects and 
    free electrons/holes as a function of electron chemical potential, mu_e, and
    temperature, T. defects is a list of instances of the Defect class, and 
    doscar is an instance of the ElectronicDOS class.
    """
    charge = 0.0
    for defect in defects:
        defect.mu_e = mu_e
        charge += (-1.)*defect.delta_e*defect.defect_concentration(T)
    charge += doscar.calc_hval(mu_e,T)
    charge -= doscar.calc_econd(mu_e,T)
    return charge

def fermi_energy(T,defects,doscar):
    """
    Calculates the electron chemical potential (fermi energy) as a function of
    temperature. 
    """
    El = doscar.Ef
    Eh = doscar.Ef
    while charge_neutrality(El,T,defects,doscar) < 0:
        El -= BOLTZCONST*T
    while charge_neutrality(Eh,T,defects,doscar) > 0:
        Eh += BOLTZCONST*T
    mu_e = brentq(charge_neutrality,El,Eh,args=(T,defects,doscar),xtol=1e-7)
    return mu_e

def main(args):
    """
    Run program logic.
    """
    electron_file = str(args[0])
    defect_file = str(args[1])
    compound_file = str(args[2])

    # Read in electronic of perfect crystal.
    e_dos = ElectronicDOS(electron_file,'ezvasp')

    # Calculate atomic chemical potentials.
    phases = input_energies(compound_file)
    chempots = ChemicalPotential(phases)
    mu = chempots.mu
    
    # Read in list defects to consider.
    defects = input_defect(defect_file)
    for defect in defects:
        defect.mu = mu

    # Calculate the electronic chemical potential at several temperatures, as
    # well as the concentrations of atomic defects, conduction band electrons,
    # and valence band holes.
    temps = np.linspace(0,2000,21)
    mu_e_array = []
    c_defects = []
    for T in temps:
        mu_e = fermi_energy(T,defects,e_dos)
        mu_e_array.append(mu_e)
        c_defects.append([])
        for defect in defects:
            defect.mu_e = mu_e
            c_defects[-1].append(defect.defect_concentration(T))
        c_defects[-1].append(e_dos.calc_econd(mu_e,T))
        c_defects[-1].append(e_dos.calc_hval(mu_e,T))

    # Print out the defect concentrations.
    for i in range(len(temps)):
        line = str(temps[i])
        line += " " + str(mu_e_array[i])
        for j in range(len(c_defects[i])):
            line += " " + str(c_defects[i][j])
        print line

if __name__ == "__main__":
    main(sys.argv[1:])
