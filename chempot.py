#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve
from math import *
from mix import *
import matplotlib.pyplot as plt
from SAFT import *


#Eq. A.5

def mu_Ass(T,dens_num,mix):

		num_c = len(mix.xcomp)

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

		#d[g_jk(djk)hs] / d[p_i]
		deriv_rdf = np.zeros((num_c,num_c,num_c))
		for i in range(0,num_c):
			for j in range(0,num_c):
				for k in range(0,num_c):
					term1 = d[i]**3 /(1-zeta[3])**2
					term2 =  d[i]**2 / (1-zeta[3])**2 + 2*(d[i]**3)*zeta[2] / (1-zeta[3])**3
					term2 = term2 * 3*d[j]*d[k]/(d[j]+d[k]) 
					term3a = 2.0 * (d[i]**2) * zeta[2] / (1-zeta[3])**3 
					term3b = 3.*(d[i]**3) * (zeta[2]**2) / (1-zeta[3])**4
					term3 = term3a + term3b
					term3 = term3 * 2.*(d[j]*d[k]/(d[j]+d[k])**2)
					deriv_rdf[j,k,i] = pi*mix.m[i]*( term1 + term2 + term3 ) / 6.0 #Eq. A.5
					print deriv_rdf[j,k,i]

		#Partial Derivative of Association Strength 
		for i in range(0,num_c):
			indx1= -1
			for j in range(0,num_c):
				for a in range(0,mix.num_assocs[i]):
					indx1 += 1
					indx2 = -1
					for k in range(0,num_c):
						for b in range(0,mix.num_assocs[k]):
							indx2 += 1
							if i == k: #same component
								kappa = kappa[i][a,b]
								epsilon = mix.eps_ass[i][a,b]
							else:
								k1 = mix.kappa[i].max()
								k2 = mix.kappa[k].max()
								e1 = mix.eps_ass[i].max()
								e2 = mix.eps_ass[k].max()
								kappa = sqrt(k1*k2) * (sqrt(mix.sigma[i]*mix.sigma[k])/(0.5*(mix.sigma[i]+mix.sigma[k])))**3
								epsilon = 0.5*(e1+e2)

							djk = (d[j] +d[k])/2.0
							deriv_delta[indx1,indx2,i] = djk**3 * deriv_rdf[j,k,i] * kappa*(exp(epsilon/T) - 1)
						












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

EtOH1 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.6)
EtOH2 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.4)
mix = Mix(EtOH1,EtOH2)
mu_Ass(T,dens_num,mix)

