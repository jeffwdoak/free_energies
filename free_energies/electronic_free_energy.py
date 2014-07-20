#!/usr/bin/python

# electronicdos.py v1.3 8-06-2013 Jeff Doak jeff.w.doak@gmail.com
import numpy as np
from scipy.optimize import brentq
import sys
import subprocess

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

    Class Methods:
    """

    def __init__(self,input_,format=None):
        if isinstance(input_,str):
            try:
                input_ = open(input_,'r')
            except IOError:
                print "Error reading input file."
                print "Program will now exit!"
                sys.exit(1)
        if isinstance(input_,file):
            if format == "ezvasp":
                self.read_ezvasp_dos(input_)
            else:
                self.read_doscar(input_)
            	#nelec = subprocess.Popen("grep NELECT OUTCAR",
            	#    shell=True,stdin=None,stdout=subprocess.PIPE).communicate()[0]
            	#self.nelec = int(float(nelec.split()[2]))
            self.get_bandgap()
        # Calculate finite temperature properties
        self.temp = np.linspace(0,2000,21)
        #self.temp[0] = 1.
        self.mu_e = np.zeros_like(self.temp)
        self.num_e = np.zeros_like(self.temp)
        self.num_h = np.zeros_like(self.temp)
        self.E_el = np.zeros_like(self.temp)
        self.S_el = np.zeros_like(self.temp)
        self.F_el = np.zeros_like(self.temp)
        # Calculate E_el_0
        self.E_el_0 = None
        tol = 1e-5
        for i in range(len(self.temp)):
            if self.temp[i] < tol:
                self.mu_e[i] = self.e_fermi
                self.E_el[i] = 0.0
                self.S_el[i] = 0.0
                self.num_e[i] = 0.0
                self.num_h[i] = 0.0
            elif self.temp[i] > tol:
                self.mu_e[i] = self.calc_mu_e(self.temp[i])
                if self.E_el_0 == None:
                    self.E_el_0 = self.calc_E_el(self.mu_e[i],self.temp[i])
                self.num_e[i] = self.n(self.mu_e[i],self.temp[i])
                self.num_h[i] = self.p(self.mu_e[i],self.temp[i])
                self.E_el[i] = (self.calc_E_el(self.mu_e[i],self.temp[i]))
                self.S_el[i] = self.calc_S_el(self.mu_e[i],self.temp[i])

        self.E_el[1:] = self.E_el[1:] - self.E_el_0
        self.F_el = self.E_el - self.temp*BOLTZCONST*self.S_el

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
	self.n_dos = int(line[2])
        self.e_fermi = float(line[3])
        energy = []; dos_tot = []; dos_spin = []
        #for line in input_:
	for i in range(self.n_dos):
            line = input_.readline().split()
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
        #self.dos_tot = np.array(dos_tot)
        #self.dos_spin = np.array(dos_spin)

    def read_ezvasp_dos(self,input_):
        """
        Reads an ezvasp-formatted dos.out file to get the electronic density of
        states. The argument input_ is assumned to be a file object.
        """
        nions = subprocess.Popen("grep NIONS OUTCAR",
                shell=True,stdin=None,stdout=subprocess.PIPE).communicate()[0]
        self.n_atoms = int(float(nions.split()[-1]))
        self.e_min = 0.0
        line = input_.readline().split()
        self.nelec = int(float(line[0]))
        self.step_size = float(line[1])
        self.scale = float(line[2])
        energy = []; dos_tot = []
        i = 0
        for line in input_:
            line = line.split()
            dos_tot.append(float(line[0]))
            energy.append(float(i)*self.step_size)
            i += 1
        self.energy = np.array(energy)
        self.dos_tot = np.array(dos_tot)
        self.dos_spin = np.zeros_like(self.dos_tot) # Change this for spin-polar
        #self.dos_spline = UnivariateSpline(self.energy,self.dos_tot)
        self.e_max = self.energy[-1]
        # Find the 0 Kelvin 'Fermi Energy' using ATAT's method
        ne = 0.0
        for i in range(len(self.dos_tot)):
            ne += self.dos_tot[i]*self.step_size
            e_fermi = self.energy[i]
            if ne >= self.nelec:
                break
        self.e_fermi = e_fermi

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

    def n(self,mu_e,T):
        """
        Calculate the intrinsic number of conduction electrons per atom at an
        electron chemical potential mu_e and temperature T.
        """
        def fermi(x,args):
            mu = args[0]; T = args[1]
            return 1./(np.exp((x-mu)/(BOLTZCONST*T))+1.)
        n = self.sum_dos(fermi,mu_e,self.e_max,args=(mu_e,T))
        return n

    def p(self,mu_e,T):
        """
        Calculate the intrinsic number of valence holes per atom at an electron
        chemical potential of mu_e and temperature T.
        """
        def fermi(x,args):
            mu = args[0]; T = args[1]
            return 1./(np.exp((mu-x)/(BOLTZCONST*T))+1.)
        p = self.sum_dos(fermi,self.e_min,mu_e,args=(mu_e,T))
        return p

    def charge_neutrality(self,mu_e,args):
        """
        Condition for charge neutrality for intrinsic doping in a perfect
        semiconductor. This function should be overwritten for a more
        complicated case. 
        """
        T = args  # Args could also include atomic chemical potentials.
        return self.p(mu_e,T) - self.n(mu_e,T)

    def bracket_mu_e(self,args):
        """
        Function to find a bracket around the root of charge_neutrality.
        Returns a tuple containing the lower and upper bounds for the root.
        """
        a = (self.vbm+self.cbm)/2.
        b = (self.vbm+self.cbm)/2.
        if self.band_gap < 1e-5:
            dx = self.step_size*5.
        else:
            dx = (self.cbm-self.vbm)/10.
        g = np.sqrt(2.)
        fa = self.charge_neutrality(a,args)
        fb = self.charge_neutrality(b,args)
        while True:
            a = a - dx
            fa = self.charge_neutrality(a,args)
            if np.sign(fa) != np.sign(fb):
                break
            b = b + dx
            fb = self.charge_neutrality(b,args)
            if np.sign(fa) != np.sign(fb):
                break
            dx = (b - a)*g
        return a,b

    def calc_mu_e(self,temp):
        """
        Calculate the electron chemical potential at temperature temp using the
        condition of charge neutrality.
        """
        lower,upper = self.bracket_mu_e(temp)
        mu_e = brentq(self.charge_neutrality,lower,upper,args=(temp))
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
        E = self.sum_dos(fermi_energy,self.e_min,self.e_max,args=(mu_e,T))
        #E_0 = self.sum_dos(energy,self.e_min,self.e_fermi,args=None)
        return E
    
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
                return -f*np.log(f)-(1.-f)*np.log(1.-f)
            else:
                return 0.0
        S = self.sum_dos(weight,self.e_min,self.e_max,args=(mu_e,T))
        return S

def main(argv):
    edos_name = str(argv[0])
    if edos_name == "dos.out":
	edos = ElectronicDOS(edos_name,"ezvasp")
    else:
        edos = ElectronicDOS(edos_name)
    print "T_(K), mu_e_(eV), n_e_(#/atom), n_h_(#/atom), E_el_(eV/atom), S_el_(eV/K/atom), F_el_(eV/atom)"
    for i in range(len(edos.temp)):
        print edos.temp[i],edos.mu_e[i],edos.num_e[i],edos.num_h[i],edos.E_el[i],edos.S_el[i],edos.F_el[i]
    sys.exit()

if __name__ == "__main__":
    main(sys.argv[1:])
