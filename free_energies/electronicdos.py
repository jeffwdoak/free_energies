#!/usr/bin/python

# electronicdos.py v0.4 1-01-2012 Jeff Doak jeff.w.doak@gmail.com
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
from scipy.optimize import fsolve
import sys

BOLTZCONST = 8.617e-5 #eV/K

class ElectronicDOS:
    """
    Class to calculate equilibrium carrier concentrations, as well as
    equilibrium thermodynamic properties of an electronic density of states.

    Class constants of ElectronicDOS:
    - BOLTZCONST - Boltzmann's Constant (eV/K)

    Instance attributes of ElectronicDOS:
    - n_atoms - number of atoms of the unit cell from which the DOS was
        calculated
    - energy - numpy array of energies at which DOS is calculated (eV)
    - dos_tot - numpy array of density of spin up and spin down allowed electron
        states at each energy in the array energy (# states/eV/atom)
    - dos_spin - numpy array of the difference in density between spin up and
        spin down states (# states/eV/atom)
    - e_min - minimum energy in numpy array energy (eV)
    - e_max - maximum energy in numpy array energy (eV)
    - e_fermi - zero Kelvin fermi energy for the electronic DOS (eV)
    - step_size - energy difference between two consecutive points in the DOSCAR
        file (eV)
    - vbm - valence band maximum, set to e_fermi for metals (eV)
    - cbm - conduction band minimum, set to e_fermi for metals (eV)
    - band_gap -  band gap around the fermi energy, zero for metals (eV)
    - temp - numpy array of temperatures at which to calculate equilibrium
        electron chemical potentials, electron concentrations, and hole
        concentrations (K)
    - mu_e - numpy array of electron chemical potentials calculated at each
        temperature in temp (eV)
    - num_e - numpy array of equilibrium electron concentrations calculated at
        each temperature in temp (# e's/atom)
    - num_h - numpy array of equilibrium hole concentrations calculated at each
        temperature in temp (# h's/atom)
    - E_el - numpy array of electronic energy calculated at each temperature
        in temp (eV/atom)
    - S_el - numpy array of electronic entropy calculated at each temperature
        in temp (kB/atom)
    - F_el - numpy array of electronic free energy calculated at each
        temperature in temp (eV/atom)
    """

    def __init__(self,input_):
        if isinstance(input_,str):
            try:
                input_ = open(input_,'r')
            except IOError:
                print "Error reading input file."
                print "Program will now exit!"
                sys.exit(1)
        if isinstance(input_,file):
            self.read_doscar(input_)
        self.get_bandgap()
        # Calculate finite temperature properties
        self.temp = np.linspace(0,2000,21)
        self.mu_e = np.zeros_like(self.temp)
        self.num_e = np.zeros_like(self.temp)
        self.num_h = np.zeros_like(self.temp)
        self.E_el = np.zeros_like(self.temp)
        self.S_el = np.zeros_like(self.temp)
        self.F_el = np.zeros_like(self.temp)
        for i in range(len(self.temp)):
            self.mu_e[i] = self.calc_mu_e(self.temp[i])
            self.num_e[i] = self.n(self.mu_e[i],self.temp[i])
            self.num_h[i] = self.p(self.mu_e[i],self.temp[i])
            self.E_el[i] = self.calc_E_el(self.mu_e[i],self.temp[i])
            self.S_el[i] = self.calc_S_el(self.mu_e[i],self.temp[i])
            self.F_el[i] = self.E_el[i] - self.temp[i]*BOLTZCONST*self.S_el[i]

    def read_doscar(self,input_):
        """
        Reads in a doscar file to grab the density of states as a function of
        energy. The argument input_ is assumed to be a file object.
        """
        self.n_atoms = int(input_.readline().split()[0])
        # Discard header information
        for i in range(4):
            input_.readline()
        # Read in Fermi Energy
        line = input_.readline().split()
        self.e_max = float(line[0])
        self.e_min = float(line[1])
        self.e_fermi = float(line[3])
        energy = []; dos_tot = []; dos_spin = []
        for line in input_:
            line = line.split()
            energy.append(float(line[0]))
            if len(line) == 3:
                dos_tot.append(float(line[1]))  # DOS includes spin up and down
                dos_spin.append(0.0)
            elif len(line) == 5:
                dos_tot.append(float(line[1])+float(line[2]))
                dos_spin.append(float(line[1])-float(line[2]))
        self.energy = np.array(energy)
        self.dos_tot = np.array(dos_tot)/float(self.n_atoms)
        self.dos_spin = np.array(dos_spin)/float(self.n_atoms)
        self.dos_spline = UnivariateSpline(self.energy,self.dos_tot,s=0)

    def get_bandgap(self):
        """
        Finds the band gap of a DOS around the fermi energy.
        """
        self.step_size = self.energy[1] - self.energy[0]
        i = 0
        not_found = True
        while not_found:
            if self.energy[i] < self.e_fermi and self.dos_tot[i] > 1e-3:
                bot = self.energy[i]
            elif self.energy[i] > self.e_fermi and self.dos_tot[i] > 1e-3:
                top = self.energy[i]
                not_found = False
            i += 1
        if top - bot < 2*self.step_size:
            self.vbm = self.cbm = self.e_fermi
            self.band_gap = 0.0
        else:
            self.vbm = bot; self.cbm = top
            self.band_gap = top - bot

    def shift_energy(self,new_ref):
        """
        Change the reference energy for all of the energy attributes.
        """
        self.energy = self.energy - new_ref
        self.e_min = self.e_min - new_ref
        self.e_max = self.e_max - new_ref
        self.e_fermi = self.e_fermi - new_ref
        self.vbm = self.vbm - new_ref
        self.cbm = self.cbm - new_ref
        self.mu_e = self.mu_e - new_ref

    def sum_dos(self,weight,start,end,args=None):
        """
        Sums the density of states, dos, in the energy range [start,end], weighted
        by the function weight, which takes as inputs energy and args.
        """
        flag = False
        sum = 0.
        for i in range(len(self.energy)):
            if flag:
                sum += self.step_size*self.dos_tot[i]*weight(
                        self.energy[i],args)
                if self.energy[i] > end:
                    break
            elif self.energy[i] >= start:
                flag = True
        return sum

    def integrate_dos(self,weight,start,end,args=None,threshold=0.1):
        """
        Takes numpy arrays containing the energy and dos and integrates them over
        the range [start,end] with the weighting function weight. Weight should take
        as an argument the integrated energy and a list of other arguements args.
        """
        def integrand(x,weight,args):
            return self.dos_spline(x)*weight(x,args)
        result = quad(
                integrand,start,end,args=(weight,args),full_output=1,limit=350)
        integral = result[0]
        error = result[1]
        #if error > integral*threshold:
        #    print "Numerical integration error is greater than"
        #    print str(threshold)+" of the integrated value."
        #    sys.exit(1)
        return integral

    def n(self,mu_e,T):
        """
        Calculate the intrinsic number of conduction electrons per atom at an
        electron chemical potential mu_e and temperature T.
        """
        def fermi(x,args):
            mu = args[0]; T = args[1]
            return 1./(np.exp((x-mu)/(BOLTZCONST*T))+1.)
        #n = self.integrate_dos(fermi,self.cbm,self.e_max,args=(mu_e,T))
        n = self.sum_dos(fermi,self.cbm,self.e_max,args=(mu_e,T))
        return n

    def p(self,mu_e,T):
        """
        Calculate the intrinsic number of valence holes per atom at an electron
        chemical potential of mu_e and temperature T.
        """
        def fermi(x,args):
            mu = args[0]; T = args[1]
            return 1./(np.exp((mu-x)/(BOLTZCONST*T))+1.)
        #p = self.integrate_dos(fermi,self.e_min,self.vbm,args=(mu_e,T))
        p = self.sum_dos(fermi,self.e_min,self.vbm,args=(mu_e,T))
        return p

    def charge_neut2(self,mu_e,args):
        T = args[0]; n_elec = args[1]
        n_sum = self.sum_dos(fermi,self,e_start,self.e_end,args=(mu_e,T))
        return n_elec - n_sum

    def charge_neutrality(self,mu_e,args):
        """
        Condition for charge neutrality for intrinsic doping in a perfect
        semiconductor. This function should be overwritten for a more
        complicated case. 
        """
        T = args  # Args could also include atomic chemical potentials.
        return self.p(mu_e,T) - self.n(mu_e,T)

    def calc_mu_e(self,temp):
        """
        Calculate the electron chemical potential at temperature temp using the
        condition of charge neutrality.
        """
        mu_e = fsolve(self.charge_neutrality,self.e_fermi,args=(temp))
        return mu_e

    def calc_E_el(self,mu_e,T):
        """
        Calculate the electronic energy at a temperature T and electron chemical
        potential mu_e.
        """
        def energy(x,args):
            return x
        def fermi_energy(x,args):
            mu = args[0]; T = args[1]
            if x-mu < -30.0*BOLTZCONST*T:
                return x
            elif x-mu > 30.0*BOLTZCONST*T:
                return 0.0
            else:
                return x/(np.exp((x-mu)/(BOLTZCONST*T))+1.)
        #E = self.integrate_dos(fermi_energy,self.e_min,self.e_max,args=(mu_e,T))
        #E_0 = self.integrate_dos(
        #        fermi_energy,self.e_min,self.e_max,args=(mu_e,T))
        E = self.sum_dos(fermi_energy,self.e_min,self.e_max,args=(mu_e,T))
        E_0 = self.sum_dos(energy,self.e_min,self.e_fermi,args=None)
        return E - E_0
    
    def calc_S_el(self,mu_e,T):
        """
        Calculate the electronic entropy at an electron chemical potential mu_e
        and temperature T.
        """
        def weight(x,args):
            mu = args[0]; T = args[1]
            x = (x - mu)/(BOLTZCONST*T)
            f = 1.0/(np.exp(x)+1)
            if f > 1e-5 and (1.0 - f) > 1e-5:
                return f*np.log(f)+(1.-f)*np.log(1.-f)
            else:
                return 0.0
            #f = -np.log(np.exp(x)+1)/(np.exp(x)+1)
            #f += -np.log(np.exp(-x)+1)/(np.exp(-x)+1)
            #return f
        #S = self.integrate_dos(weight,self.e_min,self.e_max,args=(mu_e,T))
        S = self.sum_dos(weight,self.e_min,self.e_max,args=(mu_e,T))
        return S

def fermi_dirac_dist(x,args):
    """
    Calculates the Fermi-Dirac distribution for an energy x, temperature
    args[0], and electron chemical potential args[1].
    """
    T = args[0]; mu = args[1]
    return 1./(np.exp((x-mu)/(BOLTZCONST*T))+1.)

def test2(argv):
    doscar = ElectronicDOS(open(str(argv[0]),'r'))
    T = 500
    #n = doscar.integrate_dos(
    #        fermi_dirac_dist,doscar.cbm,doscar.e_max,args=(T,doscar.e_fermi))
    p = doscar.p(doscar.e_fermi,T)
    print p

def test3(argv):
    doscar = ElectronicDOS(open(str(argv[0]),'r'))
    print doscar.temp
    print doscar.num_e
    print doscar.num_h
    print doscar.E_el
    print doscar.S_el
    print doscar.F_el

def test1(argv):
    import matplotlib.pyplot as plt
    doscar = ElectronicDOS(open(str(argv[0]),'r'))
    plt.plot(doscar.energy,doscar.dos_tot)
    plt.show()

def main(argv):
    import matplotlib.pyplot as plt
    doscar = open(str(argv[0]))
    e_fermi,energy,n_tot,n_spin = read_doscar(doscar)
    plt.plot(energy,n_tot)
    if len(argv) > 1:
        doscar2 = open(str(argv[1]))
        e_fermi2,energy2,n_tot2,n_spin2 = read_doscar(doscar2)
        plt.plot(energy2,n_tot2)
    plt.show()

if __name__ == "__main__":
    import sys
    test3(sys.argv[1:])
