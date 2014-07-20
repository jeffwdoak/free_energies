#!/software/enthought/python/epd-7.1-2-rh5-x86_64/bin/python

# bin_mix_workspace.py v0.2 5/7/2012 Jeff Doak jeff.w.doak@gmail.com

# This script uses the Phase and BinaryMixingModel classes to load up the
# thermodynamic properties of mixing for several compositions in a pseudo-binary
# system. The goal is that this can be called interactively to then perform
# specific tasks with the binary mixing model of the binary phases.

import scipy as sp
from phase import Phase
from binarymixingmodel import *
import matplotlib.pyplot as plt

# Test for creating a phase and plotting thermodynamic data about this phase.
def phase_test():
    # Read in compositon and 0-K energy.
    energyfile = "energies_0_K"
    efile = open(energyfile,'r')
    efile.readline()
    temp = efile.readline().split()
    comp = float(temp[0])
    energy = float(temp[1])
    efile.close()
    
    inputfile = "0_thermo"
    phase = Phase(comp,energy,inputfile)
    
    # Plot free energy vs temperature.
    plt.plot(phase.T,phase.F-phase.F[0])
    plt.xlabel('T (K)')
    plt.ylabel('F_vib (meV/atom)')
    plt.show()

def zero_k_energy_plot(comp,energy0):
    eform = binary_mixing_property(energy0,comp,2)
    plt.plot(comp,eform,'o')
    plt.xlabel('x_PbTe')
    plt.ylabel('Formation Energy (meV/cation)')
    plt.show()

def phases_test():
    energyfile = "energies_0_K"
    vibfiles = ["0_thermo","1.54_thermo","1.3_thermo","2.3_thermo","1_thermo"]
    efile = open(energyfile,'r')
    efile.readline()
    comp = sp.array([])
    energy0 = sp.array([])
    for line in efile:
        comp = sp.append(comp,float(line.split()[0]))
        energy0 = sp.append(energy0,float(line.split()[1]))
    efile.close()
    phase = []
    for i in range(len(comp)):
        phase.append(Phase(comp[i],energy0[i],vibfiles[i]))

    # Text output of data
    for i in range(len(phase[0].T)):
        text = str(phase[0].T[i])
        for j in range(len(phase)):
            text += " "+str(phase[j].E_thermal[i]*2.)
        print text
    sys.exit()

    # Graphical output of data
    for i in range(len(phase[0].T)):
        text = str(phase[0].T[i])
        for j in range(len(phase)):
            text += " "+str(phase[j].E_vib[i])
        text +="\n"
        print text

    for i in range(len(phase)):
        plt.plot(phase[i].T,phase[i].E_thermal*2.)
    plt.xlabel('T (K)')
    plt.ylabel('S_vib (k_B/cation)')
    plt.show()
    
def mix_prop_test():
    energyfile = "energies_0_K"
    vibfiles = ["0_thermo","1.54_thermo","1.3_thermo","2.3_thermo","1_thermo"]
    efile = open(energyfile,'r')
    efile.readline()
    comp = sp.array([])
    energy0 = sp.array([])
    for line in efile:
        comp = sp.append(comp,float(line.split()[0]))
        energy0 = sp.append(energy0,float(line.split()[1]))
    efile.close()
    phase = []
    for i in range(len(comp)):
        phase.append(Phase(comp[i],energy0[i],vibfiles[i]))
    vib_prop = []
    for i in range(len(phase)):
        vib_prop.append(phase[i].E_thermal)
    vib_mix_prop = binary_mixing_property(vib_prop,comp,2)

    # Text output of data
    for i in range(len(phase[0].T)):
        text = str(phase[0].T[i])
        for j in range(len(phase)):
            text += " "+str(vib_mix_prop[j,i])
        print text
    sys.exit()

    # Graphical output of data
    for i in range(len(phase)):
        plt.plot(phase[i].T,vib_mix_prop[i])
    plt.xlabel('T (K)')
    plt.ylabel('Delta E_vib (meV/cation)')
    plt.show()

def comp_dependence():
    energyfile = "energies_0_K"
    vibfiles = ["0_thermo","1.54_thermo","1.3_thermo","2.3_thermo","1_thermo"]
    efile = open(energyfile,'r')
    efile.readline()
    comp = sp.array([])
    energy0 = sp.array([])
    for line in efile:
        comp = sp.append(comp,float(line.split()[0]))
        energy0 = sp.append(energy0,float(line.split()[1]))
    efile.close()
    phase = []
    for i in range(len(comp)):
        phase.append(Phase(comp[i],energy0[i],vibfiles[i]))
    energy_0k = []
    vib_energy = []
    vib_entropy = []
    helm = []
    for i in range(len(phase)):
        energy_0k.append(phase[i].E_0)
        vib_energy.append(phase[i].E_vib)
        vib_entropy.append(phase[i].S_vib)
        helm.append(phase[i].F)
    mix_energy0 = binary_mixing_property(energy0,comp,2)
    mix_energy_0K = binary_mixing_property(energy_0k,comp,2)
    vib_mix_energy = binary_mixing_property(vib_energy,comp,2)
    vib_mix_entropy = binary_mixing_property(vib_entropy,comp,2)
    vib_mix_helm = binary_mixing_property(helm,comp,2)
    t=100
    print comp
    print mix_energy0
    print mix_energy_0K
    print phase[0].T[t]
    print vib_mix_energy[:,t]
    print vib_mix_entropy[:,t]*(-1.)*phase[0].T[t]
    print vib_mix_helm[:,t]

def model_fitting():
    energyfile = "energies_0_K"
    vibfiles = ["0_thermo","1.54_thermo","1.3_thermo","2.3_thermo","1_thermo"]
    efile = open(energyfile,'r')
    efile.readline()
    comp = sp.array([])
    energy0 = sp.array([])
    for line in efile:
        comp = sp.append(comp,float(line.split()[0]))
        energy0 = sp.append(energy0,float(line.split()[1]))
    efile.close()
    phase = []
    for i in range(len(comp)):
        phase.append(Phase(comp[i],energy0[i],vibfiles[i]))

    binmodel = BinaryMixingModel(phase,1,1,2,"S_vib")
    print binmodel.fit_vec




if __name__ == "__main__":
    model_fitting()
    comp_dependence()
    #phases_test()
    #mix_prop_test()
