#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve
from math import *
from mix import *
import matplotlib.pyplot as plt
from SAFT import *
from zfactor import *
from internal import *

def mu_Seg(T,dens_num,mix):
	#Fetch a0, z0, u0
	a0, epsm, a01, a02 = Helmholtz_Seg(T,dens_num,mix)
	zseg, z0 = Z_Seg(T,dens_num,mix)
	u0 = u_Seg(T,dens_num,mix)

	d = [0,0]
	num_c = 2

	#Calculate temp-dependent segment diameter for each compound
	for i in range(0,num_c):
		d[i] = HS_Diameter(T,mix.m[i],mix.sigma[i],mix.epsilon[i])


	#Reduced Density, eta
	eta_term = 0
	for i in range(0,num_c):
		eta_term += mix.xcomp[i]*mix.m[i]*d[i]**3 #Eq. 33 (d is modified)
	eta = (pi * dens_num * eta_term)/6.0 

	#Calculate HS contribution
	a_hs = (4*eta - 3*eta**2) / (1-eta)**2 #Eq. 31

	#Calculate sigma/epsilon shit
	sigma_mix = np.zeros((num_c,num_c))
	epsilon_mix = np.zeros((num_c,num_c))
	for i in range(0,num_c):
		for j in range(0,num_c):
			"""Epsilon_IJ takes a mixing factor called mix.k, which follows Peng-Robison mixing rules, 
			Chapman et al says this can reduce to 1, so I'm going to keep it here and try to figure it out later"""
			epsilon_mix[i,j] = 1.0 * sqrt(mix.epsilon[i]*mix.epsilon[j]) #Eq. 6 (but note my comment above)
			sigma_mix[i,j] = (mix.sigma[i] + mix.sigma[j]) /2.0 #Eq. 7

	denominator = 0
	for i in range(0,num_c):
		denominator += mix.xcomp[i]*mix.m[i]

	sigmaX_numerator = 0
	sigmaXepsilonX_numerator = 0
	for i in range(0,num_c):
		for j in range(0,num_c):
			sigmaX_numerator +=  mix.xcomp[i]*mix.xcomp[j]*mix.m[i]*mix.m[j]*sigma_mix[i,j]**3
			sigmaXepsilonX_numerator +=  mix.xcomp[i]*mix.xcomp[j]*mix.m[i]*mix.m[j]*epsilon_mix[i,j]*sigma_mix[i,j]**3



	sigmaXcubed = sigmaX_numerator / (denominator**2) #Eq. 4 
	sigmaX = sigmaXcubed ** 0.333333333
	sigmaXepsilonX = sigmaXepsilonX_numerator / (denominator**2) #Eq. 5
	epsilonX = sigmaXepsilonX / sigmaXcubed

	sum1 = np.zeros((num_c))
	sum2 = np.zeros((num_c))
	sum3 = np.zeros((num_c))
	term1and3 = np.zeros((num_c))
	term2 = np.zeros((num_c))
	museg = np.zeros((num_c))



	for i in range(0,num_c):
		for j in range(0,num_c):
			sum1[i] += mix.xcomp[j]*mix.m[j]*sigma_mix[i,j]**3
			sum2[i] += mix.xcomp[j]*mix.m[j]
			sum3[i] += mix.xcomp[j]*mix.m[j]*epsilon_mix[i,j]*sigma_mix[i,j]**3

		term1and3[i] = 2*sum1[i]/(sigmaXcubed*sum2[i]) 
		term2[i] = 2*sum3[i]/(epsilonX*sigmaXcubed*sum2[i])

		museg[i] = a0 + (z0-1.0)*(term1and3[i] - 1.0) + u0*(term2[i] - term1and3[i])
	return museg










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



