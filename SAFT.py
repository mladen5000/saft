#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve
from math import *
from mix import *



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
			print zeta[i]

		#Calculate the RDF
		RDF = np.zeros((num_c,num_c))
		for i in range(0,num_c):
			for j in range(0,num_c):
				#Eq. 25 and Eq. 26
				term1 = 1.0 / (1.0 - zeta[3])
				term2 = 3.*d[i]*d[j] / (d[i]+d[j])  * zeta[2] / (1. - zeta[3])**2
				term3 = 2.*( d[i]*d[j] / (d[i]+d[j]))**2  * zeta[2]**2 / (1-zeta[3])**3
				RDF[i,j] = term1 + term2 + term3
		print "THE RDF IS:", RDF

		
		#Chain contribution to Helmholtz Energy
		a_chain = 0
		for i in range(0,num_c):
				a_chain += mix.xcomp[i] * (1. - mix.m[i]) * log(RDF[i,i])

		return a_chain

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
		return a_seg

def Xa_Calculation(Xa,T,dens_num,mix,d,RDF):
		
		#Association Strength
		delta = np.zeros((num_assocs,num_assocs))
		for i in range(0,num_c):
			for j in range(0,num_assocs):
				for k in range(0,num_c):
					for l in range(0,num_assocs):

						if i==k:
							kappa = mix.kappa[i][j,l]
							epsilon = mix.epsilon[i][j,l]
						else:
							#These are some abritrary mixing rules
							#FIXMELATER 1. need to work out 'max', 2. need to understand these mixing rules
							k1 = max(mix.kappa[i])
							k2 = max(mix.kappa[k])
							e1 = max(mix.epsilon[i])
							e2 = max(mix.epsilon[k])
							kappa = sqrt(k1*k2) * (sqrt(mix.sigma[i]*mix.sigma[k])/(0.5*(mix.sigma[i]+mix.sigma[k])))**3
							epsilon = 0.5*(e1 + e2)

						dij = (d[i]+d[k])/2.0
						delta[j,l] = dij**3 * RDF[i,k] * kappa * exp( (epsilon/T) - 1.0 )


		"""
		count1 = 0
		for i in range(0,num_assocs):
			count1 = count1 + 1
			count2 = 0
			for j in range(0,num_assocs):
				kappa_current = kappa.item((i,j))
				epsilon = eps_ass.item((i,j))
				delta[i,j] = d**3 * RDF * kappa_current * (exp(epsilon/T) - 1)
		"""

		#Calculate Equation for Xa
		res = np.zeros((num_assocs,1))
		count1 = 0
		for i in range(0,num_assocs):
			count1 = count1 + 1
			count2 = 0
			sumtotal = 0

			for j in range(0,num_assocs):
				sumtotal = sumtotal + xcomp*dens_num*Xa[j]*delta[i,j]

			res[i] = Xa[i] - (1+sumtotal)**(-1)

		#Calculate the Jacobian
		J = np.zeros((num_assocs,num_assocs))

		for i in range(0,num_assocs):
			sumtotal = 0
			for j in range(0,num_assocs):
				sumtotal = sumtotal + xcomp*dens_num*Xa[j]*delta[i,j]
			for j in range(0,num_assocs):
				J[i,j] = J[i,j] + xcomp*dens_num*delta[i,j] / (1+sumtotal)**2
			J[i,i] = J[i,i] + 1

		#Only works for n = 2 right now
		output = [res.item(0),res.item(1)]
		return output


def Helmholtz_Ass(T,dens_num,mix):
		"""Provides the Association contribution to Helmholtz energy
		"""

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
			print zeta[i]

		#Calculate the RDF
		RDF = np.zeros((num_c,num_c))
		for i in range(0,num_c):
			for j in range(0,num_c):
				#Eq. 25 and Eq. 26
				term1 = 1.0 / (1.0 - zeta[3])
				term2 = 3.*d[i]*d[j] / (d[i]+d[j])  * zeta[2] / (1. - zeta[3])**2
				term3 = 2.*( d[i]*d[j] / (d[i]+d[j]))**2  * zeta[2]**2 / (1-zeta[3])**3
				RDF[i,j] = term1 + term2 + term3
		print "THE RDF IS:", RDF

		#Initial Values to guess for X^i
		#Only works for n=2 right now
		x0 = [0.91, 0.91]

		Xa = fsolve( Xa_Calculation , x0 , args=(T,dens_num,m,sigma,kappa,eps_ass,xcomp,num_assocs,d,RDF) )

		#Finally Calculate A_Assoc
		A_assoc = 0
		for i in range(0,num_assocs):
			A_assoc = A_assoc + log(Xa[i]) - Xa[i]/2 
		A_assoc = A_assoc + 0.5*num_assocs
		return A_assoc


def Sum_Helmholtz(c1):

		#General Required Input Parameters
		T=313.0 #Temperature (Kelvins)
		dens_num=0.001 #Density number (Not molar density)
		m=c1.m #Avg number of segments per chain 
		sigma=c1.sigma #Temperature Independent diameter (Angstroms)
		epsilon = c1.epsilon #LJ Interaction energy (Kelvins)

		#Association Parameters
		num_assocs= c1.num_assocs
		#Arrays [ [AA AB], [BA, BB] ]
		kappa = c1.kappa #Association Volume (Dimensionless)
		eps_ass = c1.eps_ass #Association Energy (Kelvins)

		#Composition
		xcomp=1

		#Call Necessary functions
		a_assoc = Helmholtz_Ass(T,dens_num,m,sigma,kappa,epsilon,xcomp,num_assocs,eps_ass)
		a_seg 	= Helmholtz_Seg(T,dens_num,c1)
		a_chain = Helmholtz_Chain(T,dens_num,m,sigma,epsilon,xcomp)
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
	
"""Initialize it"""
T=313.0 #Temperature (Kelvins)
dens_num=0.001 #Density number (Not molar density), given by molar density x N_Avogadro
m=2.457 #Avg number of segments per chain 
sigma=3.044 #Temperature Independent diameter (Angstroms)
epsilon = 213.48 #LJ Interaction energy (Kelvins), ((technically epsilon/k))
num_assocs=2
kappa = np.array([[0., 0.0292],[0.0292,0]]) #Association Volume (Dimensionless)
eps_ass = np.array([[0, 2619],[2619,0]]) #Association Energy (Kelvins)

EtOH1 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.5)
EtOH2 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.5)
my_mixture = Mix(EtOH1,EtOH2)
print my_mixture.kappa[0][0,1]
#print Helmholtz_Chain(T,dens_num,my_mixture)


#Sum_Helmholtz(EtOH)



		
