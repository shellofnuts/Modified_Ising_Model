# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 12:11:08 2021

@author: Robert Clampett
"""
# Change Log
"""
Updated the lattice altering to do sweeps of the lattice instead of random points. This should
still satisfy the Markov principle.

"""
"""
Removed some extraneous code that was calculating values that could be extracted after the simulation.

Trying out the random vector change. Seeing what effect changing p0_new only after acceptance of the vector 
has on computation time and on the results.
"""

# Imports

import numpy as np
import pandas as pd

# Variable Declaration

anisotropy_Exch = 0 # The magnetic exchange anisotropy
anisotropy_Axis = 1 # The energy associated with alignment to the easy axis

init_latticeSize = 40 # Defining the lattice size as nxn

init_warmupSteps = 200 # How many times to alter the random, initialised lattice 
simulating_Steps = 300# How many times to alter the lattice at each Temperature
sampling_Frequency = 1 # Save the lattice after ever n alters

starting_Temp = 1.0 # Temperatire to start the simulation at
stopping_Temp = 1.5 # Temperatire to finish
interval_Temp = 1.0005 # Temperature step

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
	# The random vector should be new for each grid point as the chance to change direction must be
	# independent for each spin, i.e. if the vector was updated after acceptance, we would have 
	x, y, z = np.shape(lattice)
	check = 0
	
	p0_new = random_Spin()
	
	for m in range(x):
		for n in range(y):
			if check == 1:
				check = 0
				p0_new = random_Spin()
			
			p0 = lattice[m][n]
			# Use of modulo to simulate boundary conditions
			neighbours = [lattice[m-1][n], lattice[m][n-1], lattice[(m+1) % x][n], lattice[m][(n+1) % y]]
	
			if accept_Change(T, p0, p0_new, neighbours):
				lattice[m][n] = p0_new
				check = 1

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
	
	test_mags = []
	
	for i in range(test_Steps):
		test_Lattice = alter_Lattice(test_Lattice, sample_Temp)
		if i % sample_Freq == 0:
			test_mags.append(magnetisation(test_Lattice))
	
	# We return a set of snapshots of the sample at the designated temperature
	# This allows for analysis to be done at a higher function
	return test_mags, test_Lattice


def MCrun(lattice_Size, warmup_Steps, test_Steps, sample_Freq, start_Temp, end_Temp, interval_Temp):
	# Here we implement the MC simulation across multiple temperatures to find
	# the critical phenomena/phase transitions
	
	global anisotropy_Axis, anisotropy_Exch
#	temp_range = np.flip(temp_range)
	
#	test_lattice = initiate_Lattice(lattice_Size)
	test_lattice = np.zeros((lattice_Size, lattice_Size, 3))
	for i in range(lattice_Size):
		for j in range(lattice_Size):
			test_lattice[i][j] = [0,0,1]
	
	test_lattice = warmup(test_lattice, warmup_Steps, start_Temp)

	magnetic_dist = pd.DataFrame()
	
	t = start_Temp
	
	while t < end_Temp:
		test_samples, test_lattice = Metropolis(test_lattice, test_Steps, sample_Freq, t)
		print(t)
		magnetic_dist[str(t)] = test_samples
		t *= interval_Temp
		
	magnetic_dist.to_csv('{}x{}_spin_distribution_TempSweep_RandVec.csv'.format(lattice_Size,lattice_Size), index=False)

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

	
# Execution	

MCrun(init_latticeSize, init_warmupSteps, simulating_Steps, sampling_Frequency, starting_Temp, stopping_Temp, interval_Temp)