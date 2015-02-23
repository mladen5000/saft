#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve
from math import *
from mix import *
import matplotlib.pyplot as plt
from SAFT import *
from zfactor import *
from internal import *
from museg import *

def mu_Total(T,dens_num,mix):
	muass = mu_Ass(T,dens_num,mix)
	muchain = mu_Chain(T,dens_num,mix)
	museg = mu_Seg(T,dens_num,mix)
	mutotal = np.zeros((num_c))

	for i in range(0,num_c):
		mutotal[i] = muass[i] + muchain[i] + museg[i]
	print mutotal
	return mutotal

def binodal(T,dens_num,mix):
	chempot = mu_Total(T,dens_num,mix)
	diff = chempot[0] - chempot[1]
	print diff
	return diff





"""Initialize it"""
T=313.0 #Temperature (Kelvins)
dens_num=0.001 #Density number (Not molar density), given by molar density x N_Avogadro
m=2.457 #Avg number of segments per chain 
sigma=3.044 #Temperature Independent diameter (Angstroms)
epsilon = 213.48 #LJ Interaction energy (Kelvins), ((technically epsilon/k))
num_assocs=2
kappa = np.array([[0.0, 0.0292],[0.0292,0.0]]) #Association Volume (Dimensionless)
eps_ass = np.array([[0.0, 2619],[2619,0.0]]) #Association Energy (Kelvins)
num_c = 2
# ----------------------------------------------------------------------------------------------
#FIXMELATER Playing with the xcomp doesn't match matlab values

EtOH1 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.50)
EtOH2 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.50)
mix = Mix(EtOH1,EtOH2)

binodal(T,dens_num,mix)

