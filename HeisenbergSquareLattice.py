# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 12:11:08 2021

@author: Robert Clampett
"""

# Imports

import numpy as np
import pandas as pd
from scipy.stats import moment

# Variable Declaration

anisotropy_Exch = 0 # The magnetic exchange anisotropy
anisotropy_Axis = 5 # The energy associated with alignment to the easy axis

init_latticeSize = 80 # Defining the lattice size as nxn

init_warmupSteps = 400000 # How many times to alter the random, initialised lattice 
simulating_Steps = 400000 # How many times to alter the lattice at each Temperature
sampling_Frequency = 1000 # Save the lattice after ever n alters

starting_Temp = 2.0 # Temperatire to start the simulation at
stopping_Temp = 0.7 # Temperatire to finish
interval_Temp = 0.01 # Temperature step

# Functions

def initiate_Lattice(n):
	# Initialise an input lattice of size nxn with random spins
	# This is done before starting "cool down" procedure
	
	lattice = np.zeros((n, n, 3))
	
	for i in range(n):
		for j in range(n):
			lattice[i][j] = random_Spin()
	
	return lattice

			
def random_Spin():
	# Returns a random unit vector
	# Using a guassian distribution to choose from achieves a probability density
	# that only depends on r. Since there is no dependence on theta and phi, 
	# this will yield a random unit vector, with equal probability to point in any direction.
	spin_Vec = [np.random.normal() for i in range(3)]
	norm_Vec = np.linalg.norm(spin_Vec)
	return spin_Vec/norm_Vec


def Hamiltonian(p0, neighbours):
	# Returns the energy associated with a lattice position p0
	# It uses the nearest neighbours p1-p4 and the normalised interaction energy
	# Uses various global variables to introduce anisotropy.
	
	global anisotropy_Exch, anisotropy_Axis
	energy = 0
	# Sum up the altered dot product
	for i in range(len(neighbours)):
		energy += -1*((1-anisotropy_Exch)*(p0[0]*neighbours[i][0] + p0[1]*neighbours[i][1]) + p0[2]*neighbours[i][2])
	
	# Add the easy axis anisotropy
	energy += -2*(anisotropy_Axis)*(p0[2]**2)
	
	return energy


def acceptance_Ratio(energy, T_kb):
	return np.exp(-1*energy/(T_kb))
	

def accept_Change(T, p0, p0_n, neighbours):
	# Here we check to see if the new spin decreases the energy or not.
	
	original_E = Hamiltonian(p0, neighbours)
	new_E = Hamiltonian(p0_n, neighbours)
	
	# Three cases need to be checked
	# 1. The energy is decreased by the change.
	# 2. The energy is increased but this is acceptable due to temp.
	# 3. The energy is increased but it is not acceptable due to temp.
	
	delta_E = new_E - original_E
	
	if delta_E <= 0 :
		#print("0\n")
		return True
	elif np.random.random() <= acceptance_Ratio(delta_E, T):
		#print("1\n")
		return True
	else:
		#print("2\n")
		return False


def alter_Lattice(lattice, T):
	# Here we want to pick a random point in the lattice and see if we can randomly change it.
	
	x, y, z = np.shape(lattice)
	m, n = np.random.randint(x), np.random.randint(y)
	
	p0 = lattice[m][n]
	# Use of modulo to simulate boundary conditions
	neighbours = [lattice[m-1][n], lattice[m][n-1], lattice[(m+1) % x][n], lattice[m][(n+1) % y]]
	
	p0_new = random_Spin()
	
	if accept_Change(T, p0, p0_new, neighbours):
		lattice[m][n] = p0_new
	return lattice


def warmup(lattice, max_steps, T):
	# We will alter the lattice at a high temperature to equilibrate the lattice
	# Our initial random lattice could be unnatural.
	for i in range(max_steps):
		lattice = alter_Lattice(lattice, T)
	return lattice

	
def Metropolis(test_Lattice, test_Steps, sample_Freq, sample_Temp):
	# Here we combine the random changing with sampling
	# We do this at a constant temperature and use our acceptance probability
	# to sample around the mean in the n-dimensional state space.
	
	test_samples = []
	
	for i in range(test_Steps):
		test_Lattice = alter_Lattice(test_Lattice, sample_Temp)
		if i % sample_Freq == 0:
			test_samples.append(np.copy(test_Lattice))
	
	# We return a set of snapshots of the sample at the designated temperature
	# This allows for analysis to be done at a higher function
	return test_samples, test_Lattice


def MCrun(lattice_Size, warmup_Steps, test_Steps, sample_Freq, start_Temp, end_Temp, interval_Temp):
	# Here we implement the MC simulation across multiple temperatures to find
	# the critical phenomena/phase transitions
	
	global anisotropy_Axis, anisotropy_Exch
	
	temp_range = np.arange(end_Temp, start_Temp, interval_Temp)
	temp_range = np.flip(temp_range)
	
	MC_avg_magnetisation = []
	binder_temps2 = []
	binder_temps4 = []
	
	test_lattice = initiate_Lattice(lattice_Size)
	test_lattice = warmup(test_lattice, warmup_Steps, start_Temp)
#	test_lattice = np.zeros((lattice_Size, lattice_Size, 3))
#	for i in range(lattice_Size):
#		for j in range(lattice_Size):
#			test_lattice[i][j] = [0,0,1]
#	
#	
#	test_lattice = warmup(test_lattice, warmup_Steps, end_Temp)
	magnetic_dist = pd.DataFrame()
	
	for i in temp_range:
		sample_mag = []
#		binderParameter_Average = []
		test_samples, test_lattice = Metropolis(test_lattice, test_Steps, sample_Freq, i)
		print(i)
		length, ii, jj, kk = np.shape(test_samples)
		for j in range(length):
			sample_mag.append(magnetisation(np.copy(test_samples[j])))
#			binderParameter_Average.append(binderParameter2(np.ndarray.flatten(test_samples[j][:][:][2])))
		MC_avg_magnetisation.append(np.average((sample_mag)))
		binder_temps2.append((binderParameter2(sample_mag)))
		binder_temps4.append((binderParameter(sample_mag)))
		magnetic_dist[str(i)] = sample_mag
	
	binderTempsDF = pd.DataFrame()
	binderTempsDF['Temps'] = temp_range
	binderTempsDF['Binder Cumulant U2'] = binder_temps2
	binderTempsDF['Binder Cumulant U4'] = binder_temps4
	
	binderTempsDF.to_csv('{}x{}_binder_cumulant.csv'.format(lattice_Size,lattice_Size), index=False)
	magnetic_dist.to_csv('{}x{}_spin_distribution.csv'.format(lattice_Size,lattice_Size), index=False)

# Characteristic Functions

def magnetisation(lattice):
	# Returns the average magnetism per site in the lattice
	normVec = np.mean(np.mean(lattice, axis=0), axis=0)
	
#	x = np.average(lattice[:][:][0])
#	y = np.average(lattice[:][:][1])
#	z = np.average(lattice[:][:][2])
	
	# Return the length of the "average" magnetic spin
	return np.linalg.norm(normVec)
#	return np.linalg.norm([x, y, z])
	
def susceptibility(lattice, T):
		return (1/T)*np.var(lattice)

def binderParameter(set_of_mag):
	return 1 - moment(set_of_mag, moment=4, axis=None)/(3*moment(set_of_mag, moment=2, axis=None)**2)
	
def binderParameter2(set_of_mag):
	return moment(set_of_mag, moment=2, axis=None)/(np.average(np.abs(set_of_mag))**2)
	
	
# Execution	

MCrun(init_latticeSize, init_warmupSteps, simulating_Steps, sampling_Frequency, starting_Temp, stopping_Temp, interval_Temp)
#MCrun(2*init_latticeSize, init_warmupSteps, simulating_Steps, sampling_Frequency, starting_Temp, stopping_Temp, interval_Temp)

#MCrun(init_latticeSize, init_warmupSteps, simulating_Steps, sampling_Frequency, starting_Temp, stopping_Temp, interval_Temp)