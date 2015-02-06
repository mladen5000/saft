#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve
from math import *
from mix import *
import matplotlib.pyplot as plt
from SAFT import *


#Eq. A.5
#d[g_jk(djk)hs] / d[p_i]

def mu_Ass(T,dens_num,mix):
	num_c = len(mix.xcomp)
	print num_c
	deriv_rdf = np.zeros((num_c,num_c,num_c))
	for i in range(0,num_c):
		for j in range(0,num_c):
			for k in range(0,num_c):
				print i,j,k
				term1 = d[i]**3 /(1-zeta[3])**2
				term2 =  d[i]**2 / (1-zeta[3])**2 + 2*(d[i]**3)*zeta[2] / (1-zeta[3])**3
				term2 = term2 * (3*d[j]*d[k]/(d[j]+d[k]) 
				term3a = 2.0 * (d[i]**2) * zeta[2] / (1-zeta[3])**3 
				term3b = 3.*(d[i]**3) * (zeta[2]**2) / (1-zeta[3])**4
				term3 = term3a + term3b
				term3 = term3 * 2.*(d[j]*d[k]/(d[j]+d[k])**2)
				deriv_rdf[j,k,i] = pi*mix.m[i]*( term1 + term2 + term3 ) / 6.0












"""Initialize it"""
T=313.0 #Temperature (Kelvins)
dens_num=0.001 #Density number (Not molar density), given by molar density x N_Avogadro
m=2.457 #Avg number of segments per chain 
sigma=3.044 #Temperature Independent diameter (Angstroms)
epsilon = 213.48 #LJ Interaction energy (Kelvins), ((technically epsilon/k))
num_assocs=0
kappa = np.array([[0.0, 0.0292],[0.0292,0.0]]) #Association Volume (Dimensionless)
eps_ass = np.array([[0.0, 2619],[2619,0.0]]) #Association Energy (Kelvins)
num_c = 2
# ----------------------------------------------------------------------------------------------

EtOH1 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.5)
EtOH2 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.5)
mix = Mix(EtOH1,EtOH2)
mu_Ass(T,dens_num,mix)

