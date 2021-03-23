# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 15:12:06 2020

@author: s1708916
"""
import numpy as np
import math
import pandas as pd
import time

# Function to randomly assign n_pharma to locations on map
# One indicates presence of pharmacy, zero no pharmacy
def InitPharmacies(data, n_pharma):
    
    # Shuffle data
    data = data.sample(frac = 1)                  
    n_vals = len(data)               # Number of sites on grid
    arr1 = np.zeros(n_vals-n_pharma)    # List of empty sites
    arr2 = np.ones(n_pharma)            # List of n_pharma pharmacies

    # Concatenate arr1 and arr2, and suffle
    arr3 = np.concatenate((arr1, arr2), axis = 0)
    np.random.shuffle(arr3)                     
    
    # Assign pharmacies to sites
    data['Pharmacy'] = arr3                       
    
    # Reset dataframe indices
    data = data.sort_index()                      

    return arr3


# Randomly select a pharmacy, a, and a trial site, b
def TrialSwitchLoc(data):
    # Sample a pharmacy
    map_pharma = data[data['Pharmacy']>=1.0] 
    a = map_pharma.sample(1)
    
    # Sample a random site on grid
    b = data.sample(1)
    
    return a, b


# Find energies of all pharmacies in initial setup
def InitEnergies(data, max_dist, mean_pop, r_c, r_r, sf):
    
    # Initialise energy column in dataframe
    data['Energy'] = np.zeros(len(data))
    
    # Find indices of sites with pharmacies
    map_pharma = data[data['Pharmacy']>=1.0]
    pharma_indices = list(map_pharma.index)
    j = 0
    
    print("Start energy initialisation")
    
    # Calculate energy of all pharmacies
    for i in pharma_indices: 
        j+= 1
        x_pos = data['x'].loc[i] # x coordinate of pharmacy
        y_pos = data['y'].loc[i] # y coordinate of pharmacy

        # Calculate distance to all points on map
        data_sphere = DistSquared(data, i, x_pos, y_pos, max_dist)

        # Calculate attractive contributions of all points in pharmacy catchment area
        attractive = attraction(data_sphere, mean_pop, r_c)

        # Remove pharmacy at original site to avoid it interacting with itself
        data_sphere['Pharmacy'].loc[i] -= 1

        # Calculate repulsive contributions of all pharmacies within r_r
        repulsive = repulsion(data_sphere, r_r, mean_pop, sf)

        # Calculate energy of pharmacy
        data['Energy'].loc[i] = (attractive + 0.5 * repulsive.sum(axis = 0))

        # Move back pharmacy
        data_sphere['Pharmacy'].loc[i] += 1

        print("Initial energy calculated", j, "out of", len(pharma_indices))
    
    return data


# Calculate the distance squared of a pharmacy to surrounding data points within a radius of max_dist
def DistSquared(data, i, xpos, ypos, max_dist):
    
    # Initialise squared distance to point of interest in dataframe
    data['Dist sq'] = np.zeros(len(data))

    # Calculate distance only for data points within max_dist in x and y directions
    data_sphere = data[(data['x'] > (xpos - max_dist)) & (data['x'] < (xpos + max_dist))]
    data_sphere = data_sphere[(data_sphere['y'] > (ypos - max_dist)) & (data_sphere['y'] < (ypos + max_dist))]
    
    # Calculate the squared distance to all data points
    data_sphere['Dist sq'] = (data_sphere['x'] - xpos)**2 + (data_sphere['y'] - ypos)**2
    
    # Keep only values within a radius of max_dist of pharmacy to define the sphere of influence
    data_sphere = data_sphere[data_sphere['Dist sq'] <= max_dist**2]
    
    return data_sphere


# Find repulsive contributions of all pharmacies in sphere of influence
# Modelled by decreasing linear function with x-intercept r_r, scaled by sf
def repulsion(data_sphere, r_r, mean_pop, sf):

    # Dataframe of data points within a radius of r_r from pharmacy site
    repulsive_pharma = data_sphere[data_sphere['Dist sq'] <= r_r**2]

    # Calculate repulsive contribution of all pharmacies in repulsive_pharma
    repulsive = sf * repulsive_pharma['Pharmacy'] * ((r_r - np.sqrt(repulsive_pharma['Dist sq']))/r_r)
    
    return repulsive


# Find attractive contributions of population at all data points in sphere of influence
# Can be chosen as Fermi-Dirac or step-wise function.
def attraction(data_sphere, mean_pop, r_c):
    
    # For step function
    #step_left = data_sphere[data_sphere['Dist sq'] <= r_c**2]
    #attractive = step_left['Population'].sum(axis = 0) / mean_pop
    
    # For Fermi function
    attractive = (data_sphere['Population']/ mean_pop) * 1/(1+math.e**((np.sqrt(data_sphere['Dist sq'])-r_c)/r_c))
    attractive = attractive.sum(axis = 0)
    
    return -1 * attractive


# Adjust energies of surrounding pharmacies if trial move is accepted
def UpdateSurroundings(data, index, repulsive, plusminus):

    # Find indices of pharamacies whose energy will change
    repulsive_indices = list(repulsive.index)
    
    # Increase energies of pharmacies at new site to account for increased repulsion 
    if plusminus > 0 :
        for x in repulsive_indices :
            pharma_val = data['Pharmacy'].loc[x]
            if pharma_val >= 1 :
                data['Energy'].loc[x] += (0.5*repulsive.loc[x])/pharma_val
                
            
        
    # Reduce energies of pharmacies at original site to account for reduced repulsion 
    if plusminus < 0 :
        for x in repulsive_indices :
            pharma_val = data['Pharmacy'].loc[x]
            if pharma_val >= 1 :
                data['Energy'].loc[x] -= (0.5*repulsive.loc[x])/pharma_val
                  
    return data


# Calculate energy of pharmacy based on attractive and repulsive interactions
def pharma_E(attractive, repulsive):
    pharma_E = (attractive + 0.5 * repulsive.sum(axis=0))
    
    return pharma_E


# Calculate energy difference of pharmacy at initial and trial location 
def delta_E(attractive, repulsive, attractive_orig, repulsive_orig) :
    
    delta_E = (attractive + repulsive.sum(axis = 0)) - (attractive_orig + repulsive_orig.sum(axis = 0)) 
    
    return delta_E 
 

# Calculate system energy as the sum of energies of all pharmacies
def System_E(data):
    system_energy = (data['Energy']*data['Pharmacy']).sum(axis=0)
    
    return system_energy

# Test if system has converged
def TestConvergence(system_energy_vals, t, T, n_pharma, convergence):
    if t > n_pharma * 100 :
        # Select last 10 system energy values
        a = np.array(system_energy_vals)
        a = a[-10:]

        # If system energy variation is less than T, the system has converged
        if np.max(a) - np.min(a) <= 100 :
            convergence = True
    
    return convergence

# Sample the average occupation of each site on the grid
def AverageOccupation(data, t, n_sample, samples):
    # Sample average occupation every 100 timesteps
    data['Average Occupation'] = (data['Average Occupation']*n_sample + data['Pharmacy'] ) / (n_sample + 1)
        
    return data # or running list of values?

# Keep track of the time length of occupation
def Timescale(data, t, a_index):
    # Update average occupation time
    data['Occupation time'].loc[a_index] = (data['Occupation time'].loc[a_index] * (t-1) + (t - data['Reset time'].loc[a_index]) ) / t
    
    # Update reset time
    data['Reset time'].loc[a_index] = t
    
    return data
