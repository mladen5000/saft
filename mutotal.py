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
	print "MUASS",muass
	muchain = mu_Chain(T,dens_num,mix)
	print "MUCHAIN",muchain
	museg = mu_Seg(T,dens_num,mix)
	print "MUSEG",museg
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

def MUvT(T,dens_num,mix):
	tf = 900
	ti = 100
	arraylength = tf - ti
	mu = np.zeros(arraylength)
	x = np.zeros(arraylength)
	index = -1
	for i in range(ti,tf):
		index += 1
		temp = Z_TOTAL(i,dens_num,mix)
		mu[index] = temp
		x[index] = index
	plt.plot(x,mu)
	plt.show()





"""Initialize it"""
T=313.0 #Temperature (Kelvins)
dens_num=0.001 #Density number (Not molar density), given by molar density x N_Avogadro

"""ETOH"""
m=2.457 #Avg number of segments per chain 
sigma=3.044 #Temperature Independent diameter (Angstroms)
epsilon = 213.48 #LJ Interaction energy (Kelvins), ((technically epsilon/k))
num_assocs=2
kappa = np.array([[0.0, 0.0292],[0.0292,0.0]]) #Association Volume (Dimensionless)
eps_ass = np.array([[0.0, 2619],[2619,0.0]]) #Association Energy (Kelvins)
num_c = 2

"""CO2"""
m2 = 1.417
sigma2 = 3.172
epsilon2 = 216.08 
num_assocs2 = 0
# ----------------------------------------------------------------------------------------------
#FIXMELATER Playing with the xcomp doesn't match matlab values

EtOH1 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.50)
CO2 = Compound(sigma2,epsilon2,m2,num_assocs2,kappa,eps_ass,.50)
#EtOH2 = Compound(sigma,epsilon,m,num_assocs,kappa,eps_ass,.50)
mix = Mix(CO2,EtOH1)
MUvT(T,dens_num,mix)
#binodal(T,dens_num,mix)

