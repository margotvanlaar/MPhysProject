# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 14:48:21 2020

@author: s1708916
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dynamicsV4
import data_prep
import visualisation
import time


def main():

    
    # Read in trial grid (20 x 20 datapoints, resoultion 1km^2, pop = 1 at all points)
    pop_map = pd.read_csv('T = 100 Final Data.csv', header = 0, index_col = 0, dtype = float)
    
    # Create arrays to monitor convergence behaviour
    pop_map['Average Occupation'] = pop_map['Pharmacy']
    
    # Tracking loc Accra centre
    tracking_loc_Accra = pop_map['Pharmacy'].idxmax()
    
    # Tracking location countryside
    single_pharma_map = pop_map[pop_map['Pharmacy'] == 1]
    tracking_loc_countryside = single_pharma_map['Pharmacy'].idxmin()
    
    average_occupation_Accra = []
    average_occupation_countryside = []
    
    # Get occupation number on both sites
    value_Accra = pop_map['Average Occupation'].loc[tracking_loc_Accra]
    value_countryside = pop_map['Average Occupation'].loc[tracking_loc_countryside]
    
    # Append average occupation
    average_occupation_Accra.append(value_Accra)
    average_occupation_countryside.append(value_countryside)
          
    # Track number of samples
    samples = []
    n_sample = 1
    samples.append(n_sample)
    
    # Append initial system energy
    system_E = dynamicsV4.System_E(pop_map)
    
    # System pa#rameters
    n_sweeps = 40000    # Number of sweeps for Monte Carlo simulation
    max_dist = 31       # Maximum travel distance to pharmacy
    r_c = 15            # Cut off radius for attraction
    r_r = 30            # Cut off radius for repulsion
    mean_pop = pop_map['Population'].mean()  #Mean population for each datapoint
    T = 100             # Temperature of system
    sf = 1.1            # Scale factor
    
    #Create empty lists for energy plot
    system_energy_vals = []
    timestep = []
    

    # Start MC loop
    
    print("Start MC Loop")
    
    # Number of successful switches to be monitored
    success = 0
    
    # Condition for convergence
    convergence = False
    
    
    for t in range(n_sweeps):

        # Select random pharmacy, a, and random trial site, b
        a, b = dynamicsV4.TrialSwitchLoc(pop_map)
        a_index = a.index
        b_index = b.index
        
        # Find horizontal and vertical positions of trial pharmacy location
        bx_pos = b.iloc[0,3]
        by_pos = b.iloc[0,4]
        
        # Find horizontal and vertical positions of original pharmacy location
        ax_pos = a.iloc[0,3]
        ay_pos = a.iloc[0,4]
        
        # Remove pharmacy from its original location (ensure doesn't interact with ghost of itself)
        pop_map['Pharmacy'].loc[a_index] -= 1.0
        
        # Calculate distance squared to points in sphere of influence at trial site
        data_sphere = dynamicsV4.DistSquared(pop_map, b_index, bx_pos, by_pos, max_dist)
        
        # Calculate attractive and repulsive terms to all points in sphere of influence at trial site
        attractive = dynamicsV4.attraction(data_sphere, mean_pop, r_c)           
        repulsive = dynamicsV4.repulsion(data_sphere, r_r, mean_pop, sf)
        
        # Calculate distance squared to points in sphere of influence at original site
        data_sphere_orig = dynamicsV4.DistSquared(pop_map, a_index, ax_pos, ay_pos, max_dist)
        
        # Calculate attractive and repulsive terms to all points in sphere of influence at original site
        repulsive_orig = dynamicsV4.repulsion(data_sphere_orig, r_r, mean_pop, sf)
        attractive_orig = dynamicsV4.attraction(data_sphere_orig, mean_pop, r_c)
        
        
        # Calculate energy difference of switch
        E_diff = (attractive + repulsive.sum(axis = 0)) - (attractive_orig + repulsive_orig.sum(axis = 0)) 
        
        # Accept switch if energy change negative
        if E_diff <= 0 :
            
            # Move pharmacy to trial site
            pop_map['Pharmacy'].loc[b_index] += 1.0
            
            # Update energy of surrounding pharmacies
            pop_map = dynamicsV4.UpdateSurroundings(pop_map, b_index, repulsive, 1.0)
            pop_map = dynamicsV4.UpdateSurroundings(pop_map, a_index, repulsive_orig, -1.0)
            
            # Assign new energy to pharmacy at trial site
            pop_map['Energy'].loc[b_index] = (attractive + 0.5 * repulsive.sum(axis=0))
        
            # If no pharmacy left at original site, reduce it's energy to zero
            val = float(pop_map['Pharmacy'].loc[a_index])
            if val == 0 :
                pop_map['Energy'].loc[a_index] = 0 #-= (attractive_orig +  0.5 * repulsive_orig.sum(axis=0))
            
            
            # Calculate system energy from dataframe
            system_E = dynamicsV4.System_E(pop_map)
            
            success += 1
            
        # Accept change with Boltzmann probability if energy change positive 
        elif E_diff > 0 :
            p = math.exp(-E_diff/T)
            if np.random.rand(1,1) <= p :
                # Accept switch
                
                # Move pharmacy to trial site
                pop_map['Pharmacy'].loc[b_index] += 1.0
                
                # Update energy of surrounding pharmacies
                pop_map = dynamicsV4.UpdateSurroundings(pop_map, b_index, repulsive, 1.0)
                pop_map = dynamicsV4.UpdateSurroundings(pop_map, a_index, repulsive_orig, -1.0)
                
                # Assign new energy to pharmacy at trial site
                pop_map['Energy'].loc[b_index] = (attractive + 0.5 * repulsive.sum(axis=0))
            
                # If no pharmacy left at original site, reduce it's energy to zero
                val = float(pop_map['Pharmacy'].loc[a_index])
                if val == 0 :
                    pop_map['Energy'].loc[a_index] = 0 #-= (attractive_orig +  0.5 * repulsive_orig.sum(axis=0))
                
                
                # Calculate system energy from dataframe
                system_E = dynamicsV4.System_E(pop_map)
                success += 1       
        
            else:
                pop_map['Pharmacy'].loc[a_index] += 1.0
        
        
        # Calculate average occupation of each site
        if t % 100 == 0: 
            n_sample += 1
            samples.append(n_sample)
            pop_map = dynamicsV4.AverageOccupation(pop_map, t, n_sample, samples)
            
            value_Accra = float(pop_map['Average Occupation'].loc[tracking_loc_Accra])
            value_countryside = float(pop_map['Average Occupation'].loc[tracking_loc_countryside])
            
            average_occupation_Accra.append(value_Accra)
            average_occupation_countryside.append(value_countryside)
        
        # Total energy of system
        system_energy_vals.append(system_E)
        timestep.append(t)
        
        # Save system energy
        visualisation.Plot_SystemE(t, timestep, system_energy_vals, pop_map)
        
    # Save text files
    np.savetxt("AverageOccupationAccra.txt", average_occupation_Accra, delimiter = ",")
    np.savetxt("AverageOccupationCountryside.txt", average_occupation_countryside, delimiter = ",")
    np.savetxt("Converged system energy.txt", system_energy_vals, delimiter = ",")
    
    # Output csv file
    pop_map.to_csv("Congercence test data.csv")
main()