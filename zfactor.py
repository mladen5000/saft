#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve
from math import *
from mix import *
import matplotlib.pyplot as plt
from SAFT import *
from chempot import *



def Z_Chain(T,dens_num,mix):

		#Calculate temp-dependent segment diameter for each compound
		d=[0,0] #FIXMELATER
		for i in range(0,num_c):
			d[i] = HS_Diameter(T,mix.m[i],mix.sigma[i],mix.epsilon[i])

		#Calculate the Zeta terms for RDF 
		zeta = [0,0,0,0] 
		for i in range(0,4):
			for j in range(0,num_c):
				zeta[i] +=  mix.xcomp[j]*mix.m[j]*d[j]**i #Eq. 27
			zeta[i] *= dens_num*pi/6.0

		#Calculate the RDF
		RDF = np.zeros((num_c,num_c))
		for i in range(0,num_c):
			for j in range(0,num_c):
				#Eq. 25 and Eq. 26
				term1 = 1.0 / (1.0 - zeta[3])
				term2 = 3.*d[i]*d[j] / (d[i]+d[j])  * zeta[2] / (1. - zeta[3])**2
				term3 = 2.*( d[i]*d[j] / (d[i]+d[j]))**2  * zeta[2]**2 / (1-zeta[3])**3
				RDF[i,j] = term1 + term2 + term3

		#Determine p * d(lnRDF)/dp
		derivRDF = np.zeros(num_c)
		zchain = 0
		for i in range(0,num_c):
			term1 = zeta[3]/(1.0-zeta[3])**2 
			term2 = 1.5*d[i]*zeta[2]/(1.0-zeta[3])**2
			term3 = 3.0*d[i]*zeta[2]*zeta[3]/(1.0-zeta[3])**3
			term4 = (d[i]**2) * (zeta[2]**2) / (1.0-zeta[3])**3
			term5 = 1.5*(d[i]**2)*(zeta[2]**2)*zeta[3] / (1.0 - zeta[3])**4
			factor = 1.0/RDF[i,i]
			derivRDF[i] = factor * (term1 + term2 + term3 + term4 + term5) #Eq. A.12

		#Determine Z_chain
		for i in range(0,num_c):
			zchain += mix.xcomp[i]*(1-mix.m[i])*derivRDF[i] #Eq. A.11
		print "ZCHAIN",zchain
		return zchain

def Z_Seg(T,dens_num,mix):
		#Calculate temp-dependent segment diameter for each compound
		d=[0,0] #FIXMELATER
		dontcare, epsilon_X, a01, a02 = Helmholtz_Seg(T,dens_num,mix)
		for i in range(0,num_c):
			d[i] = HS_Diameter(T,mix.m[i],mix.sigma[i],mix.epsilon[i])


		#Reduced Density, eta
		eta_term = 0
		for i in range(0,num_c):
			eta_term += mix.xcomp[i]*mix.m[i]*d[i]**3 #Eq. 33 (d is modified)
		eta = (pi * dens_num * eta_term)/6.0 

		#Calculate Zo for Hard Spheres
		Zo_HS = (1.0 + eta + eta**2 + eta**3) / (1.0 - eta)**3 #Eq. A.15

		#Calculate reduced Density, pr
		pr = eta*(6.0/(pi*sqrt(2)))

		#Calculate Zo for Dispersion
		zo1_disp = pr * (-8.595 - 2.*(4.5424*pr) - 3.*(2.1268*pr*pr) + 4.*(10.285*pr**3) ) #Eq. A.17
		zo2_disp = pr * (-1.9075 + 2.*(9.9724*pr) - 3.*(22.216*pr*pr) + 4.*(15.904*pr**3) ) #Eq. A.18
		Tr = T/epsilon_X
		Zo_disp = zo1_disp/Tr + zo2_disp/(Tr**2) #Eq. A.16

		#Calculate Zo for the segment
		Zo_seg = Zo_HS + Zo_disp #Eq. A.14

		#Calculate Z_Segment
		sum = 0
		for i in range(0,num_c):
			sum += mix.xcomp[i]*mix.m[i]
		zseg = 1.0 + (Zo_seg - 1.)*sum #Eq. A.13
		print "ZSEG",zseg
		return zseg, Zo_seg


def Z_Assoc(T,dens_num,mix):
	num_c = 2 #FIXMELATER
	mu_ass = mu_Ass(T,dens_num,mix)
	A_ass,xa = Helmholtz_Ass(T,dens_num,mix)
	sum = 0
	for i in range(0,num_c):
		sum += mix.xcomp[i] * mu_ass[i]
	Z_ass = sum - A_ass
	print "ZASS",Z_ass
	return Z_ass

def Z_TOTAL(T,dens_num,mix):
	z_ass = Z_Assoc(T,dens_num,mix)
	z_chain = Z_Chain(T,dens_num,mix)
	z_seg,z0 = Z_Seg(T,dens_num,mix)
	Z = z_ass + z_chain + z_seg
	print Z
	return Z
	


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
Z_TOTAL(T,dens_num,mix)


