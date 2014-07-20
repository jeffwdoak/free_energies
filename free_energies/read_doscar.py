#!/usr/bin/python
import numpy as np

BOLTZCONST = 8.617e-5 #eV/K

def read_doscar(input_):
    """
    Reads in a doscar file to grab the density of states as a function of
    energy. The argument input_ is assumed to be a file object.
    """
    n_atoms = input_.readline().split()[0]
    # Discard header information
    for i in range(4):
        input_.readline()
    # Read in Fermi Energy
    line = input_.readline().split()
    e_max = float(line[0]); e_min = float(line[1]); e_fermi = float(line[3])
    energy = []; n_tot = []; n_spin = []
    for line in input_:
        line = line.split()
        energy.append(float(line[0]))
        if len(line) == 3:
            n_tot.append(float(line[1]))  # DOS includes spin up and down
            n_spin.append(0.0)
        elif len(line) == 5:
            n_tot.append(float(line[1])+float(line[2]))
            n_spin.append(float(line[1])-float(line[2]))
    energy = np.array(energy)
    n_tot = np.array(n_tot); n_spin = np.array(n_spin)
    energy = energy - e_fermi
    return e_fermi,energy,n_tot,n_spin

def get_bandgap(energy,dos):
    """
    Finds the band gap of a DOS around the fermi energy.
    """
    step_size = energy[1] - energy[0]
    i = 0
    not_found = True
    while not_found:
        if energy[i] < 0 and dos[i] > 1e-3:
            bot = energy[i]
        elif energy[i] > 0 and dos[i] > 1e-3:
            top = energy[i]
            not_found = False
        i += 1
    if top - bot < 2*step_size:
        top = bot = 0
    return bot,top

def integrate_dos(energy,dos,weight,start,end,args=None):
    """
    Takes numpy arrays containing the energy and dos and integrates them over
    the range [start,end] with the weighting function weight. Weight should take
    as an argument the integrated energy and a list of other arguements args.
    """
    from scipy.interpolate import UnivariateSpline
    from scipy.integrate import quad
    dos = UnivariateSpline(energy,dos,s=0)
    #def integrand(x,dos,weight,args):
    #    return dos(x)*weight(x,args)
    #result = quad(integrand,start,end,args=(dos,weight,args))
    result = quad(dos,start,end,limit=100)
    return result

def sum_dos(energy,dos,weight,start,end,args=None):
    """
    Sums the weighted density of states
    """
    condition = (energy - start == 0) or (energy - end == 0)
    cond_arr = np.extract(condition,energy)
    e_start = cond_arr[0]; e_end = cond_arr[1]
    sum = 0.
    for i in range(e_start,e_end):
        sum += dos[i]*weight(energy[i],args)
    return sum

def fermi_dirac_dist(x,args):
    """
    Calculates the Fermi-Dirac distribution for an energy x, temperature
    args[0], and electron chemical potential args[1].
    """
    T = args[0]; mu = args[1]
    return 1./(np.exp((x-mu)/BOLTZCONST*T)+1.)

def calc_np(argv):
    doscar = open(str(argv[0]))
    if len(argv) > 1:
        T = float(argv[1])
    else:
        T = 300.
    e_fermi,energy,n_tot,n_spin = read_doscar(doscar)
    vbm,cbm = get_bandgap(energy,n_tot)
    start = energy[0]
    end = energy[len(energy)-1]
    #n = integrate_dos(energy,n_tot,fermi_dirac_dist,cbm,end,args=(T,e_fermi))
    #p = integrate_dos(energy,n_tot,fermi_dirac_dist,start,vbm,args=(T,e_fermi))
    #p = integrate_dos(energy,n_tot,fermi_dirac_dist,start,end,args=(T,e_fermi))
    p = sum_dos(energy,n_tot,fermi_dirac_dist,start,end,args=(T,e_fermi))
    return p

def test_integrand(argv):
    doscar = open(str(argv[0]))
    if len(argv) > 1:
        T = float(argv[1])
    else:
        T = 300.
    e_fermi,energy,n_tot,n_spin = read_doscar(doscar)
    vbm,cbm = get_bandgap(energy,n_tot)
    start = energy[0]
    end = energy[len(energy)-1]
    print dos
    fermi_dirac_dist()

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
    n,p = calc_np(sys.argv[1:])
    print n,p
