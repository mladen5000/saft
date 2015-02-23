#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve
from math import *
from mix import *
import matplotlib.pyplot as plt



def HS_Diameter(T,m,sigma,epsilon):
		"""
		Calculates the Hard Sphere Diameter 
		T = Temp \K
		m = number of segments in chain
		"""

		fm = 0.0010477 + 0.025337*(m-1)/m 	#Eq. 3 of reference
		f = (1 + 0.2977*(T/epsilon)) / (1 + 0.33163*(T/epsilon) + fm*(T/epsilon)**2) 	#Eq. 2 of reference
		d = sigma*f 	#Eq. 1 of reference
		return d

def Helmholtz_Chain(T,dens_num,mix):
		"""
		Calculates Chain contribution to Helmholtz energy
		T = Temperature \K
		dens_num = Number Density \Number of molecules/Angstroms^3
		dens_num = Na * rho
		mix = Object that contains m, sigma, and epsilon for each item in the mixture
		"""
		
		num_c = 2 #FIXMELATER
		#To do this for a mixture, i'll have to loop through each item in the mixture, 
		#ie m(i), sigma(i), epsilon (i), im gonna keep it like this (1 mixture item) for now

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

		
		#Chain contribution to Helmholtz Energy
		a_chain = 0
		for i in range(0,num_c):
				a_chain += mix.xcomp[i] * (1. - mix.m[i]) * log(RDF[i,i])

		return a_chain, RDF

def Helmholtz_Seg(T,dens_num,mix):
		"""Provides the Segment contribution to Helmholtz energy
		H_SEG = H_HS + H_DISP
		"""
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

		
		#Calculate dispersive contribution
		p_r = eta*(6.0/(pi*sqrt(2)))
		a01_disp = p_r*(-8.5959 - 4.5424*p_r - 2.1268*p_r**2 + 10.285*p_r**3) #%Eq. 35 of reference
		a02_disp = p_r*(-1.9075 + 9.9724*p_r - 22.216*p_r**2 + 15.904*p_r**3) #%Eq. 36 of reference

		a_disp = (epsilonX/T) * ( a01_disp + a02_disp / (T/epsilonX) )	
		
		a_seg = a_hs + a_disp
		return a_seg, epsilonX, a01_disp, a02_disp

def Xa_Calculation(Xa,T,dens_num,mix,d,RDF):
		#Total number of associating sites
		num_assocs = mix.num_assocs[0]+mix.num_assocs[1] #FIXMELATER only for binary mixtures atm
		num_c = 2 #FIXMELATER
		
		#Association Strength
		delta = np.zeros((num_assocs,num_assocs))
		indx1 =-1 
		for i in range(0,num_c):
			for j in range(0,mix.num_assocs[0]): #FIXMELATER binary case only
				indx1 += 1
				indx2 =-1 
				for k in range(0,num_c):
					for l in range(0,mix.num_assocs[1]): #FIXMELATER binary case only
						indx2 += 1

						if i==k:
							kappa = mix.kappa[i][j,l]
							epsilon = mix.eps_ass[i][j,l]
						else:
							#These are some abritrary mixing rules
							#FIXMELATER 1. need to work out 'max', 2. need to understand these mixing rules
							k1 = mix.kappa[i].max()
							k2 = mix.kappa[k].max()
							e1 = mix.eps_ass[i].max()
							e2 = mix.eps_ass[k].max()
							kappa = sqrt(k1*k2) * (sqrt(mix.sigma[i]*mix.sigma[k])/(0.5*(mix.sigma[i]+mix.sigma[k])))**3
							epsilon = 0.5*(e1 + e2)

						dij = (d[i]+d[k])/2.0
						delta[indx1,indx2] = dij**3 * RDF[i,k] * kappa * (exp(epsilon/T) - 1.0 )
		#FIXMELATER The loops of ijkl need to fix j and l because num_assocs is an array, should be in mix to contain each

		#FIXMELATER fixed at 4
		res = np.zeros((4,1))
		indx1 = -1
		for i in range(0,num_c):
			for j in range(0,mix.num_assocs[0]):
				indx1 += 1
				indx2 = -1
				RHS_Xa = 0
				for k in range(0,num_c):
					for l in range(0,mix.num_assocs[1]):
						indx2 += 1

						RHS_Xa += mix.xcomp[k]*dens_num*Xa[indx2]*delta[indx1,indx2]


				res[indx1] = Xa[indx1] - (1. + RHS_Xa)**(-1) #Minimize the difference between Xa and RHS of Xa, should min to 0...ideally


		#FIXMELATER
		#Only works for n = fixed right now
		output = [res.item(0),res.item(1),res.item(2),res.item(3)]
		return output


def Helmholtz_Ass(T,dens_num,mix):
		"""Provides the Association contribution to Helmholtz energy
		"""
		num_c =2 

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

		#Initial Values to guess for X^i
		#Only works for num_assocs=4 right now
		x0 = [0.3,0.3,0.3,0.3]

		Xa = fsolve( Xa_Calculation , x0 , args=(T,dens_num,mix,d,RDF) )

		#Finally Calculate A_Assoc
		A_assoc = 0
		indx1 = -1
		for i in range(0,2):
			sum1 = 0
			for j in mix.num_assocs:
				indx1 += 1
				sum1 += log(Xa[indx1]) - Xa[indx1]/2.0
			A_assoc += mix.xcomp[i]*(sum1 + 0.5*mix.num_assocs[i])
		return A_assoc,Xa


def Sum_Helmholtz(T,dens_num,mix):

		#Call Necessary functions
		a_assoc,xa = Helmholtz_Ass(T,dens_num,mix)
		a_seg,epsx,a01,a02 	= Helmholtz_Seg(T,dens_num,mix)
		a_chain,RDF = Helmholtz_Chain(T,dens_num,mix)
		A = a_assoc + a_seg + a_chain
		"""
		#Print it out
		print " "
		print " -------Helmholtz-------"
		print " "
		print "A_association is =", a_assoc
		print "A_segment is =", a_seg
		print "A_chain is =", a_chain
		print " "
		print "Your Helmholtz energy is",  A
		print " "
		"""
		return A
	
"""Initialize it"""
"""
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

EtOH1 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.2)
EtOH2 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.8)
mix = Mix(EtOH1,EtOH2)


y = np.zeros(700)

T = np.arange(150,850)
i=0
for tval in T:
	y[i] = Sum_Helmholtz(tval,dens_num,mix)
	i += 1
plt.plot(T,y)
plt.show()
"""



