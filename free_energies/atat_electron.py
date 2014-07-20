
# port of atat felec.c++ code to python. Uses DOSCAR instead of dos.out
# electronicdos.py v0.5 5-18-2012 Jeff Doak jeff.w.doak@gmail.com
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
from scipy.optimize import fsolve
import sys,subprocess
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
        #self.get_bandgap()
        #self.shift_energy(self.e_min)
        # Total # electrons in calculation
        #nelec = subprocess.Popen("grep NELECT OUTCAR",
        #        shell=True,stdin=None,stdout=subprocess.PIPE).communicate()[0]
        #nelec = int(float(nelec.split()[2]))
        # ATAT ported code here:
        #ne = 0.0
        #i = 0
        #while ne < nelec and i < len(self.dos_tot):
        #    ne += self.dos_tot[i]
        #    i += 1
        #Ef = self.step_size*(i-1)+self.e_min
        #print self.e_min,self.vbm,self.cbm,self.band_gap
        #print "VASP fermi energy:",str(self.e_fermi)
        #print "ATAT fermi energy:",str(Ef)


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
        self.nelec = float(line[0])
        self.step_size = float(line[1])
        scale = float(line[2])
        energy = []; dos_tot = []
        i = 0
        for line in input_:
            line = line.split()
            dos_tot.append(float(line[0]))
            energy.append(float(i)*self.step_size)
            i += 1
        self.energy = np.array(energy)
        self.dos_tot = np.array(dos_tot)

    def calc_nelec(mu,T):
        ne=0.0
        for i in range(len(self.dos_tot)):
            ne += self.dos_tot[i]/(np.exp((self.energy[i]-mu)/(BOLTZCONST*T))+1.)
        return ne

    def calc_ncond(mu,T):
        ne = 0.0
        for i in range(len(self.dos_tot)):
            if self.energy[i] > mu:
                ne += self.dos_tot[i]/(np.exp((self.energy[i]-mu)/(BOLTZCONST*T))+1.)
        return ne

    def calc_hval(mu,T):
        nh = 0.0
        for i in range(len(self.dos_tot)):
            if self.energy[i] < mu:
                nh += self.dos_tot[i]/(np.exp((mu-self.energy[i])/(BOLTZCONST*T))+1.)
        return nh

    def calc_Selec(mu,T):
        S = 0.0
        for i in range(len(self.dos_tot)):
            f = 1.0/(np.exp((self.energy[i]-mu)/(BOLTZCONST*T))+1.)
            if f > 1e-5 and (1.0-f) > 1e-5:
                S += self.dos_tot[i]*(f*np.log(f)+(1.-f)*np.log(1.-f))
        return S

    def calc_Eelec(mu,T):
        E = 0.0
        for i in range(len(self.dos_tot)):
            e = self.energy[i]*self.step_size
            if (e-mu) < -30.*BOLTZCONST*T:
                f = 1.0
            elif (e-mu) > 30.*BOLTZCONST*T:
                f = 0.0
            else:
                f = 1./(np.exp((e-mu)/(BOLTZCONST*T))+1.)
            E += self.dos_tot[i]*e*f
        return E





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
        try:
            self.vbm = self.vbm - new_ref
            self.cbm = self.cbm - new_ref
            self.mu_e = self.mu_e - new_ref
        except:
            pass

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

def test2(argv):
    import matplotlib.pyplot as plt
    doscar1 = ElectronicDOS(open(str(argv[0]),'r'))
    for i in range(len(doscar1.dos_tot)):
        if doscar1.dos_tot[i] > 1e-2:
            e = doscar1.energy[i]
            break
    doscar1.shift_energy(e)
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
    doscar = ElectronicDOS(open(str(argv[0]),'r'))
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
    test2(sys.argv[1:])
