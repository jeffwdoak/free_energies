#!/usr/bin/python

# port of atat felec.c++ code to python. Uses DOSCAR instead of dos.out
# electronicdos.py v0.7 5-23-2012 Jeff Doak jeff.w.doak@gmail.com
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
from scipy.optimize import brentq
import sys,subprocess
#BOLTZCONST = 8.617e-5 #eV/K
BOLTZCONST = 1.380658e-23/1.60217733e-19;

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
                nelec = subprocess.Popen("grep NELECT OUTCAR",
                    shell=True,stdin=None,stdout=subprocess.PIPE).communicate()[0]
                self.nelec = int(float(nelec.split()[2]))
                self.get_bandgap()
                self.dos_tot = self.dos_tot*self.step_size
                e = self.e_min
                for i in range(len(self.dos_tot)):
                    if self.dos_tot[i] > 1e-2:
                        e = self.energy[i]
                        break
                self.shift_energy(e)
        # ATAT ported code here:
        self.atat_main()

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
        self.dos_tot = np.array(dos_tot)
        self.dos_spin = np.array(dos_spin)
        self.dos_spline = UnivariateSpline(self.energy,self.dos_tot,s=0)

    def read_ezvasp_dos(self,input_):
        """
        Reads an ezvasp-formatted dos.out file to get the electronic density of
        states. The argument input_ is a ssumed to be a file object.
        """
        #self.n_atoms - get from  poscar/outcar/contcar
        self.e_min = 0.0
        line = input_.readline().split()
        self.nelec = int(float(line[0]))
        self.step_size = float(line[1])
        self.scale = float(line[2])
        energy = []; dos_tot = []
        i = 0
        for line in input_:
            line = line.split()
            dos_tot.append(float(line[0])*self.step_size)
            energy.append(float(i)*self.step_size)
            i += 1
        self.energy = np.array(energy)
        self.dos_tot = np.array(dos_tot)

    def calc_nelec(self,mu,w):
        ne=0.0
        for i in range(len(self.dos_tot)):
            ne += self.dos_tot[i]/(np.exp((float(i)-mu)/w)+1.)
        return ne

    def calc_econd(self,mu,temp):
        mu = mu/self.step_size
        w = BOLTZCONST*temp/self.step_size
        ne = 0.0
        for i in range(len(self.dos_tot)):
            if float(i) > mu:
                ne += self.dos_tot[i]/(np.exp((float(i)-mu)/w)+1.)
        return ne

    def calc_hval(self,mu,temp):
        mu = mu/self.step_size
        w = BOLTZCONST*temp/self.step_size
        nh = 0.0
        for i in range(len(self.dos_tot)):
            if float(i) < mu:
                nh += self.dos_tot[i]/(np.exp((mu-float(i))/w)+1.)
        return nh

    def calc_Selec(self,mu,w):
        S = 0.0
        for i in range(len(self.dos_tot)):
            f = 1.0/(np.exp((float(i)-mu)/w)+1.)
            if (f > 1e-5 and (1.0-f) > 1e-5):
                S += -self.dos_tot[i]*(f*np.log(f)+(1.-f)*np.log(1.-f))
        return S

    def calc_Eelec(self,mu,w,dE):
        E = 0.0
        for i in range(len(self.dos_tot)):
            e = dE*float(i)
            if (e-mu) < -30.*w:
                f = 1.0
            elif (e-mu) > 30.*w:
                f = 0.0
            else:
                f = 1./(np.exp((e-mu)/w)+1.)
            E += self.dos_tot[i]*e*f
        return E

    def atat_charge_neut(self,dE,temp):
        """
        Old ATAT version calculation of chemical potential:
        """
        ne_tol = 1e-7*self.nelec
        El = self.Ef
        Eh = self.Ef
        while self.calc_nelec(El/dE,BOLTZCONST*temp/dE) > self.nelec:
            El += -BOLTZCONST*temp
        while self.calc_nelec(Eh/dE,BOLTZCONST*temp/dE) < self.nelec:
            Eh += BOLTZCONST*temp
        while True:
            mu = (Eh+El)/2.
            ne = self.calc_nelec(mu/dE,BOLTZCONST*temp/dE)
            if ne < self.nelec:
                El = mu
            else:
                Eh = mu
            if not (np.abs(ne-self.nelec) > ne_tol):
                break
        ne = self.calc_econd(mu,temp)
        nh = self.calc_hval(mu,temp)
        return mu,ne,nh

    def charge_neutrality(self,dE,temp):
        """
        Calculate electron chemical potentials based on charge neutrality
        between the number of conduction electrons and valence holes.
        """
        ne_tol = 1e-7*self.nelec
        El = self.Ef
        Eh = self.Ef
        while self.calc_econd(El,temp) > self.calc_hval(El,temp):
            El += -BOLTZCONST*temp
        while self.calc_econd(Eh,temp) < self.calc_hval(Eh,temp):
            Eh += BOLTZCONST*temp
        while True:
            mu = (Eh+El)/2.
            ne = self.calc_econd(mu,temp)
            nh = self.calc_hval(mu,temp)
            if ne < nh:
                El = mu
            else:
                Eh = mu
            if not (np.abs(ne-nh) > ne_tol):
                break
        return mu,ne,nh

    def charge_neutrality2(self,dE,temp):
        """
        Calculate electron chemical potentials based on charge neutrality
        between the number of conduction electrons and valence holes.
        """
        ne_tol = 1e-7*self.nelec
        El = self.Ef
        Eh = self.Ef
        while self.calc_hval(El,temp) - self.calc_econd(El,temp) < 0:
            El += -BOLTZCONST*temp
        while self.calc_hval(Eh,temp) - self.calc_econd(Eh,temp) > 0:
            Eh += BOLTZCONST*temp
        while True:
            mu = (Eh+El)/2.
            ne = self.calc_econd(mu,temp)
            nh = self.calc_hval(mu,temp)
            if nh - ne > 0:
                El = mu
            else:
                Eh = mu
            if not (np.abs(nh-ne) > ne_tol):
                break
        return mu,ne,nh

    def fsolve_neutrality(self,dE,temp,mu0):
        """
        Calculate electron chemical potentials based on charge neutrality
        between electrons and holes, using fsolve to find mu_e.
        """
        ne_tol = 1e-7
        def neutrality(mu,self,temp):
            ne = self.calc_econd(mu,temp)
            nh = self.calc_hval(mu,temp)
            return (ne - nh)
        # Find bounds for brentq
        El = self.Ef
        Eh = self.Ef
        while self.calc_econd(El,temp) > self.calc_hval(El,temp):
            El += -BOLTZCONST*temp
        while self.calc_econd(Eh,temp) < self.calc_hval(Eh,temp):
            Eh += BOLTZCONST*temp
        # Run optimization using brentq
        mu = brentq(neutrality,El,Eh,args=(self,temp),xtol=ne_tol,maxiter=100)
        ne = self.calc_econd(mu,temp)
        nh = self.calc_hval(mu,temp)
        return mu,ne,nh

    def atat_main(self):
        ne = 0.0
        dE = float(self.step_size)
        i = 0
        while ne<self.nelec and i < len(self.dos_tot):
            ne += self.dos_tot[i]
            i += 1
        self.Ef = dE*(i-1.) + self.e_min

        ne_tol = 1e-7*self.nelec
        E0 = None
        zero_tolerance = 1e-3

        temps = np.linspace(0,2000,21)
        mu_e = np.zeros_like(temps)
        E_el = np.zeros_like(temps)
        S_el = np.zeros_like(temps)
        F_el = np.zeros_like(temps)
        num_e = np.zeros_like(temps)
        num_h = np.zeros_like(temps)

        for i in range(len(temps)):
            if temps[i] < zero_tolerance:
                mu_e[i] = self.Ef
                E_el[i] = 0.
                S_el[i] = 0.
                F_el[i] = 0.
                num_e[i] = 0.
                num_h[i] = 0.
            else:
                #mu,ne,nh = self.atat_charge_neut(dE,temps[i])
                mu,ne,nh = self.charge_neutrality2(dE,temps[i])
                #mu,ne,nh = self.fsolve_neutrality(dE,temps[i],mu_e[i-1])
                mu_e[i] = mu
                num_e[i] = ne
                num_h[i] = nh
                S_el[i] = self.calc_Selec(mu/dE,BOLTZCONST*temps[i]/dE)
                if temps[i] > 0.0:
                    if E0 == None:
                        E0 = self.calc_Eelec(mu,BOLTZCONST*(temps[1]-temps[0]),dE)
                    E_el[i] = self.calc_Eelec(mu,BOLTZCONST*temps[i],dE) - E0
                F_el[i] = E_el[i] - BOLTZCONST*temps[i]*S_el[i]
        self.temps = temps
        self.mu_e = mu_e
        self.E_el = E_el
        self.S_el = S_el
        self.F_el = F_el
        self.num_e = num_e
        self.num_h = num_h

    def get_bandgap(self):
        """
        Finds the band gap of a DOS around the fermi energy.
        """
        self.step_size = self.energy[1] - self.energy[0]
        tol = 1e-3
        i = 0
        not_found = True
        while not_found:
            if self.energy[i] < self.e_fermi and self.dos_tot[i] > tol:
                bot = self.energy[i]
            elif self.energy[i] > self.e_fermi and self.dos_tot[i] > tol:
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
        try:
            self.mu_e = self.mu_e - new_ref
        except:
            pass

def test2(argv):
    import matplotlib.pyplot as plt
    doscar1 = ElectronicDOS(open(str(argv[0]),'r'))
    doscar2 = ElectronicDOS(open(str(argv[1]),'r'),format="ezvasp")
    plt.plot(doscar1.energy,doscar1.dos_tot)
    plt.plot(doscar2.energy,doscar2.dos_tot)
    plt.show()

def test3(argv):
    doscar = ElectronicDOS(open(str(argv[0]),'r'))
    print doscar.temp
    print doscar.num_e
    print doscar.num_h
    print doscar.E_el
    print doscar.S_el
    print doscar.F_el

def test1(argv):
    if len(argv) < 2:
        doscar = ElectronicDOS(open(str(argv[0]),'r'))
    else:
        doscar = ElectronicDOS(open(str(argv[0]),'r'),format=str(argv[1]))
    for i in range(len(doscar.temps)):
        print doscar.temps[i],doscar.mu_e[i],doscar.E_el[i],doscar.S_el[i],doscar.F_el[i],doscar.num_e[i],doscar.num_h[i]
    #plt.plot(doscar.energy,doscar.dos_tot)
    #plt.show()

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
    test1(sys.argv[1:])
