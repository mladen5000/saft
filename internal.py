#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve
from math import *
from mix import *
import matplotlib.pyplot as plt
from SAFT import *
from chempot import *
from zfactor import *



def u_Seg(T,dens_num,mix):
	"""
	Calculate the internal energy
	"""

	#Fetch epsilon_mix, a01_disp, a02_disp
	blah,epsilon_m, a01, a02 = Helmholtz_Seg(T,dens_num,mix)

	#Calculate weighted avg of segment number
	m_prom = 0;
	for i in range(0,num_c):
		m_prom += mix.m[i]*mix.xcomp[i]

	#Calculate f(m)
	fm = 0.0010477 + 0.025337*(m_prom-1)/m_prom

	#Calculate (T/d)*(dD/dT)
	Tr = T/epsilon_m
	term1 = (1.0 - fm*(Tr**2) ) / (1 + .33163*Tr + fm*(Tr**2) )   
	term2 = 1.0/ ( 1 + 0.2977*Tr)
	rhs_term = term1 - term2 #Eq. A21

	#Calculate the first deriv from Eq. A20
	lhs_term = a01/Tr + 2*a02/(Tr**2)

	#Fetch Z0_Segment
	zseg, z0seg = Z_Seg(T,dens_num,mix)

	#Calculate u0
	u0 = lhs_term - 3.*(z0seg - 1.)*rhs_term #Eq. A19
	return u0


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

EtOH1 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.5)
EtOH2 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.5)
mix = Mix(EtOH1,EtOH2)

