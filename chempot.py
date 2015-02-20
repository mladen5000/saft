#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve
from math import *
from mix import *
import matplotlib.pyplot as plt
from SAFT import *



def mu_Ass(T,dens_num,mix):

		num_c = len(mix.xcomp)
		tot_num_assocs = mix.num_assocs[0] + mix.num_assocs[1]#FIXMELATER for nonbinary

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

		#Calculate the RDF term
		chain, RDF = Helmholtz_Chain(T,dens_num,mix)

		#d[g_jk(djk)hs] / d[p_i]
		#Eq. A.5
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
					term3 *= 2.*( (d[j]*d[k])/(d[j]+d[k]))**2 
					deriv_rdf[j,k,i] = pi*mix.m[i]*( term1 + term2 + term3 ) / 6.0 

		#Partial Derivative of Association Strength 
		#Eq. A.4
		deriv_delta = np.zeros((tot_num_assocs,tot_num_assocs,num_c))
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
								kappa = mix.kappa[i][a,b]
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


							chain, RDF = Helmholtz_Chain(T,dens_num,mix)



		#dXa/dpi Calculation
		dXa, Xa = dXa_Calculation(T,dens_num,mix,deriv_delta,deriv_rdf,d,RDF)

		#Convert dXa into a matrix
		dXa_matrix = np.zeros(( sum(mix.num_assocs),num_c ))
		for i in range(0,num_c):
			indx1 = -1;
			for j in range(0,num_c):
				for a in range(0,num_c):
					indx1 += 1
					dXa_matrix[indx1,i] = dXa[(i-1) * sum(mix.num_assocs) + indx1]

		#Put it all together for Mu_Association
		term1 = np.zeros((num_c))
		indx1 = -1
		for i in range(0,num_c):
			sum1 = 0
			for a in range(0,mix.num_assocs[i]):
				indx1 += 1
				sum1 += log(Xa[indx1]) - Xa[indx1]/2.0

			term1[i] = sum1 + 0.5*mix.num_assocs[i]

		term2 = np.zeros((num_c))
		for i in range(0,num_c):
			indx1 = -1
			sum1 = 0
			for j in range(0,num_c):
				for a in range(0,mix.num_assocs[j]):
					indx1 += 1
					sum1 += dens_num*mix.xcomp[j]* (dXa_matrix[indx1,i]*(1.0/Xa[indx1]-0.5))

			term2[i] = sum1

		muass = np.zeros((num_c))

		for i in range(0,num_c):
			muass[i] = term1[i] + term2[i]
		print muass



def dXa_Calculation(T,dens_num,mix,deriv_delta,deriv_rdf,d,RDF):
		""" This is a set of equations used by fsolve to find the available root
		"""
		A_assoc, Xa = Helmholtz_Ass(T,dens_num,mix)
		num_c = 2 #FIXMELATER
		A = np.zeros( ( num_c*sum(mix.num_assocs),num_c*sum(mix.num_assocs) ) )
		B= np.zeros( num_c*sum(mix.num_assocs) )
		dimen1 = sum(mix.num_assocs)
		delta1 = np.zeros(( dimen1 , dimen1 ))
		indx3 = -1
		for i in range(0,num_c):
			indx1 = -1
			for j in range(0,num_c):
				for a in range(0,mix.num_assocs[j]):
					indx1 += 1
					indx3 += 1
					indx2 = -1
					sum1 = 0
					for k in range(0,num_c):
						for b in range(0,mix.num_assocs[k]):
							indx2 += 1



							if j == k:
								kappa = mix.kappa[i][a,b]
								epsilon = mix.eps_ass[i][a,b]
							else:
								k1 = mix.kappa[i].max()
								k2 = mix.kappa[k].max()
								e1 = mix.eps_ass[i].max()
								e2 = mix.eps_ass[k].max()
								kappa = sqrt(k1*k2) * (sqrt(mix.sigma[i]*mix.sigma[k])/(0.5*(mix.sigma[i]+mix.sigma[k])))**3
								epsilon = 0.5*(e1+e2)


							#CALCULATE DELTA1

							dij = (d[i]+d[k])/2.0
							delta1[indx1,indx2] = dij**3 * RDF[i,k] * kappa * (exp(epsilon/T) - 1.0 )

							sum1 += dens_num*mix.xcomp[k]*(Xa[indx2]*deriv_delta[indx1,indx2,i])
							s1 =indx1 + (i-1) * sum(mix.num_assocs) #for example, skip every 4
							s2 =indx2 + (i-1) * sum(mix.num_assocs) 

							A[s1,s2] += Xa[indx1]**2 * dens_num * mix.xcomp[k] * delta1[indx1,indx2]

					sum2 = 0
					for k in range(0,mix.num_assocs[i]):
						if i == j:
							kappa = mix.kappa[j][a,k]
							epsilon = mix.eps_ass[j][a,k]
						else:
							k1 = mix.kappa[j].max()
							k2 = mix.kappa[i].max()
							e1 = mix.eps_ass[j].max()
							e2 = mix.eps_ass[i].max()
							kappa = sqrt(k1*k2) * (sqrt(mix.sigma[i]*mix.sigma[k])/(0.5*(mix.sigma[i]+mix.sigma[k])))**3
							epsilon = 0.5*(e1+e2)

						#Double Check the 2 deltas, dont think theyre the same
						dij = ( d[i] + d[j] )/2.0
						delta = dij**3 * RDF[j,i] * kappa * (exp(epsilon/T) - 1.0)
						x =sum(mix.num_assocs[0:i]) + k
						sum2 += Xa[ x ] * delta

					A[indx3,indx3] += 1.0
					B[indx3] = -(Xa[indx1])**2*(sum1 + sum2)

		#B/A' in python form
		invA= np.linalg.inv(A.T)
		res = np.dot(B,invA)

		return res,Xa
		

					
















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
mu_Ass(T,dens_num,mix)

