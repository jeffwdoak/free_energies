#!/software/enthought/python/epd-7.1-2-rh5-x86_64/bin/python
#quaternary-phases.py v1.0 4-16-2012 Jeff Doak jeff.w.doak@gmail.com

#This program calculates tie lines on an isothermal section (T=0K) of a 
#ternary phase diagram for Na-Pb-Te.
import sys
import scipy as sp
from scipy import linalg
from scipy.optimize import fmin_slsqp
import matplotlib as mp

def gibbs(phase_vec,energies,comp_vec,stoich):
    """Calculates the Gibbs Free Energy of a system consisting of
a number of phases each having a fraction in the vector phase_vec.
The energies of each phase are given in the vector energies.
Currently assumes T=0K, P=0kBar, no external fields."""
    gibbs = sp.dot(energies,phase_vec)
    return gibbs

def mass_cons(phase_vec,energies,comp_vec,stoich):
    """Function to calculate the composition of a system with
phase fractions phase_vec and phase stoichiometries stoich and
compare it with a specified composition comp_vec. Returns the
difference between calculated and desired composition vectors."""
    comp_prime = sp.dot(stoich,phase_vec)
    difference = comp_vec - comp_prime
    return difference

def min_gibbs(phase_guess,energies,comp_vec,stoich):
    """Minimizes the total Gibbs Free Energy of the system
at a composition comp_vec with respect to the fraction of each
phase present. Returns the phase fractions that minimize the
Gibbs Free Energy."""
    args = (energies,comp_vec,stoich)
    btuple = (0.0,1.0)
    bounds = []
    for i in range(len(phase_guess)):
        bounds.append(btuple)
    iprint = 0
    phase_vec = fmin_slsqp(gibbs,phase_guess,f_eqcons=mass_cons,bounds=bounds,args=args,iprint=iprint)
    return phase_vec

#Functions for plotting tie lines and relative energies.
def coords(comp_vec):
    """Converts Gibbs triangle coordinates to cartesian
coordinates for plotting."""
    cart_vec = sp.array((comp_vec[0]/2. + comp_vec[1],comp_vec[0]*sp.sqrt(3)/2))
    return cart_vec

#num_elements = 3
#Read in compound data
#Data format: Name x_1 x_2 x_3 E_formation
if len(sys.argv) > 1:
    datafile = sys.argv[1]
else:
    datafile = "data"
compounds = []
file = open(datafile,"r")
for line in file:
    line = line.split()
    temp = [line[0]]
    for i in range(1,len(line)):
        temp.append(float(line[i]))
    compounds.append(temp)
file.close()
num_compounds = len(compounds)
num_elements = len(compounds[0])-2

#Create energies vector and stoichiometry matrix
energies = sp.zeros(num_compounds)
stoich = sp.zeros((num_elements,num_compounds))
for i in range(num_compounds):
    energies[i] = compounds[i][-1]
    stoich[:,i] = sp.array(compounds[i][1:1+num_elements])
#Loop over every pair of compounds and determine if tie-line exists
#between them
tol = 1E-6
text = "{"
for i in range(num_compounds):
    for j in range(i+1,num_compounds):
        for k in range(j+1,num_compounds):
            for l in range(k+1,num_compounds):
                phase_vec = sp.zeros(num_compounds)
                phase_vec[i] = 0.25
                phase_vec[j] = 0.25
                phase_vec[k] = 0.25
                phase_vec[l] = 0.25
                comp_vec = sp.dot(stoich,phase_vec)
                phase_min = min_gibbs(phase_vec,energies,comp_vec,stoich)
                phase_diff = linalg.norm(phase_vec-phase_min)
                #print phase_diff
                if phase_diff < tol:
                    stoich2 = sp.zeros((4,4))
                    stoich2[0,:] = stoich[:,i]
                    stoich2[1,:] = stoich[:,j]
                    stoich2[2,:] = stoich[:,k]
                    stoich2[3,:] = stoich[:,l]
                    energy2 = sp.zeros(4)
                    energy2[0] = energies[i]; energy2[1] = energies[j]
                    energy2[2] = energies[k]; energy2[3] = energies[l]
                    mu = sp.dot(sp.linalg.inv(stoich2),energy2)
                    print compounds[i][0]+"-"+compounds[j][0]+"-"+compounds[k][0]+"-"+compounds[l][0]
                    print "mu_Na= "+str(mu[0])+" mu_Pb= "+str(mu[1])+" mu_S= "+str(mu[2])+" mu_Te= "+str(mu[3])
                    print

#Plot tie lines on a Gibbs triangle
#Triangle has vertices: (0,0),(1,0),(1/2,Sqrt(3)/2)
