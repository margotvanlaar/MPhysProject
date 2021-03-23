# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 14:48:21 2020

@author: s1708916
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dynamicsV4 as dynamics
import data_prep
import visualisation


def main():
   
    # Open GIS population map data
    pop_map = pd.read_excel('gha_pd.xlsx', header = 0, dtype = float)
    
    # Change coordinate system to centred Cartesian
    data_prep.DegreeToDistance(pop_map)
    data_prep.CentredDistance(pop_map)
    
    # Lower resolution of data
    pop_map = data_prep.CoarseData(pop_map, 5)
    
    # Calculate mean population at each data point
    mean_pop = pop_map['Population'].mean()
    
    # System parameters
    n_pharma = 2600              # Number of pharmacies
    n_sweeps = n_pharma * 500    # Number of sweeps for Monte Carlo simulation
    r_c = 15                     # Cut off radius for attraction
    r_r = 30                     # Cut off radius for repulsion
    max_dist = r_r + 5                # Maximum travel distance to pharmacy
    T = 100                      # Temperature of system
    sf = 1.1 
    
    #Create empty lists for energy plot
    system_energy_vals = []
    timestep = []
    
    # Set population to 1 at all points
    #pop_map['Population'] = np.ones(len(pop_map))
    
    # Randomly allocate specified number of pharmacies to a location
    pop_map = dynamics.InitPharmacies(n_pharma)
    
    # Initialise energies of all pharmacies   
    pop_map = dynamics.InitEnergies(pop_map, max_dist, mean_pop, r_c, r_r, sf)
    
    # Append initial system energy
    system_E = dynamics.System_E(pop_map)
    
    
    # Start Monte Carlo loop
    
    print("Start MC Loop")
    
    # Number of successful switches to be monitored
    success = 0
    
    # Condition for convergence
    convergence = False
    
    
    for t in range(n_sweeps):

        # Select random pharmacy, a, and random trial site, b
        a, b = dynamics.TrialSwitchLoc(pop_map)
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
        
        # Calculate distance squared to data points in sphere of influence (catchment area) at trial site
        data_sphere = dynamics.DistSquared(pop_map, b_index, bx_pos, by_pos, max_dist)
        
        # Calculate attractive and repulsive terms to all points in sphere of influence at trial site
        attractive = dynamics.attraction(data_sphere, mean_pop, r_c)           
        repulsive = dynamics.repulsion(data_sphere, r_r, mean_pop, sf)

        # Calculate distance squared to data points in sphere of influence at original site
        data_sphere_orig = dynamics.DistSquared(pop_map, a_index, ax_pos, ay_pos, max_dist)
        
        # Calculate attractive and repulsive terms to all points in sphere of influence at original site
        attractive_orig = dynamics.attraction(data_sphere_orig, mean_pop, r_c)
        repulsive_orig = dynamics.repulsion(data_sphere_orig, r_r, mean_pop, sf)
        
        
        # Calculate energy difference of switch
        E_diff = dynamics.delta_E(attractive, repulsive, attractive_orig, repulsive_orig)

        # Accept switch if energy change negative
        if E_diff <= 0 :
            
            # Move pharmacy to trial site
            pop_map['Pharmacy'].loc[b_index] += 1.0
            
            # Update energy of surrounding pharmacies
            pop_map = dynamics.UpdateSurroundings(pop_map, b_index, repulsive, 1.0)
            pop_map = dynamics.UpdateSurroundings(pop_map, a_index, repulsive_orig, -1.0)
            
            # Assign new energy to pharmacy at trial site
            pop_map['Energy'].loc[b_index] = dynamics.pharma_E(attractive, repulsive)
        
            # If no pharmacy left at original site, reduce it's energy to zero
            val = float(pop_map['Pharmacy'].loc[a_index])
            if val == 0 :
                pop_map['Energy'].loc[a_index] = 0
            
            
            # Calculate system energy from dataframe
            system_E = dynamics.System_E(pop_map)
            
            # Add successful switch
            success += 1
        
        # Accept change with Boltzmann probability if energy change positive 
        elif E_diff > 0 :
            
            # Calculate Boltzmann probability
            p = math.exp(-E_diff/T)

            # Accept switch with Boltzmann probability
            if np.random.rand(1,1) <= p :
                
                # Move pharmacy to trial site
                pop_map['Pharmacy'].loc[b_index] += 1.0
                
                # Update energy of surrounding pharmacies
                pop_map = dynamics.UpdateSurroundings(pop_map, b_index, repulsive, 1.0)
                pop_map = dynamics.UpdateSurroundings(pop_map, a_index, repulsive_orig, -1.0)
                
                # Assign new energy to pharmacy at trial site
                pop_map['Energy'].loc[b_index] = dynamics.pharma_E(attractive, repulsive)
            
                # If no pharmacy left at original site, reduce it's energy to zero
                val = float(pop_map['Pharmacy'].loc[a_index])
                if val == 0 :
                    pop_map['Energy'].loc[a_index] = 0
                
                # Calculate system energy from dataframe
                system_E = dynamics.System_E(pop_map)
                
                success += 1       
        
            else:
                pop_map['Pharmacy'].loc[a_index] += 1.0
        
        
        
        # Total energy of system
        system_energy_vals.append(system_E)
        timestep.append(t)
        
        # Plot system energy
        visualisation.Plot_SystemE(t, timestep, system_energy_vals, pop_map)
        
        # Track successful switches
        #visualisation.TrackSwitches(t, success, convergence)
        
        # Test convergence of system
        convergence = dynamics.TestConvergence(system_energy_vals, t, T,n_pharma, convergence)
        if convergence == True:
            break

      
    # Plot total energy of system
    plt.plot(timestep, system_energy_vals)
    plt.title('Total system energy vs timestep')
    plt.xlabel('Timestep')
    plt.ylabel('System energy')
    plt.savefig("System energy vs timestep")
    plt.show()
        
    
    # Plot population map
    fig,ax = plt.subplots(1,1)    
    x_pop, y_pop, pop, area_pop = visualisation.population_map(pop_map)
    ax.scatter(x_pop, y_pop, s = area_pop, c = pop, cmap = plt.get_cmap("jet"), alpha=0.4)
    #plt.colorbar()
    #ax.set(xlim=[-100,700], ylim=[-50,750])
    
    # Plot pharmacy map
    x1, y1, pharma1, area1 = visualisation.pharma_map(pop_map)
    ax.scatter(x1, y1, s = area1, c = 'red', marker = 'X')
    
    # Set figure labels
    plt.xlabel("Horizontal distance (m)")
    plt.ylabel("Vertical distance (m)")
    plt.title("Ghana population map")
    plt.plot()
    fig.savefig("Ghana population plot", dpi = 1000)
    
    # Save final dataset
    pop_map.to_csv('Final.csv') 
    
    
main()
