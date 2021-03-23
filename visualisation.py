# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 09:42:55 2020

@author: s1708916
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Data for plotting population distribution
def population_map(data):
    x = data['x']
    y = data['y']
    pop = data['Population']
    area = data['Population']/1000
    return x, y, pop, area

# Data for plotting pharmacy distribution
def pharma_map(data):
    map_pharma = data[data['Pharmacy']>=1.0]
    x = map_pharma['x']
    y = map_pharma['y']
    pharma = map_pharma['Pharmacy']
    area = data['Population'].mean()/1000
    return x, y, pharma, area

# Capture frames for video
def CaptureFrames(t, i, data, label) :
    """
    if t > 160000 :
        if t % 500 == 0 :
            #pop_map.to_csv(label)
            # Plot population map
            fig,ax = plt.subplots(1,1)    
            #x_pop, y_pop, pop, area_pop = visualisation.population_map(pop_map)
            #ax.scatter(x_pop, y_pop, s = area_pop, c = pop, cmap = plt.get_cmap("terrain"), alpha=0.4)
           
            
            x1, y1, pharma1, area1 = visualisation.pharma_map(pop_map)
            ax.scatter(x1, y1, s = area1, c = pharma1, cmap = plt.get_cmap("jet") , marker = 'X')
            ax.set(xlim=[-100,700], ylim=[-50,750])
            label = "Frame {fnumber}".format(fnumber = i)
            fig.savefig(label, dpi = 1000)
            
            i += 1
    """
    if t > 160000 : # Start plotting
        if t % 500 == 0 : # Take snapshot every 500 timesteps
            fig,ax = plt.subplots(1,1)    
            
            # Plot population map
            #x_pop, y_pop, pop, area_pop = visualisation.population_map(pop_map)
            #ax.scatter(x_pop, y_pop, s = area_pop, c = pop, cmap = plt.get_cmap("terrain"), alpha=0.4)
           
            # Plot pharmacies
            x1, y1, pharma1, area1 = pharma_map(data)
            ax.scatter(x1, y1, s = area1, c = pharma1, cmap = plt.get_cmap("jet") , marker = 'X')
            ax.set(xlim=[-100,700], ylim=[-50,750])
            label = "Frame {fnumber}".format(fnumber = i)
            fig.savefig(label, dpi = 1000)
            
            i += 1 # Update label
    
    return i

# Track fraction of successful switches
def TrackSwitches(t, success, convergence):
    if t > 1000 :
        if t % 1000 == 0. :
            frac = success/(t+1)
            print("fraction of successful switches = ", frac )
        # Convergence condition
            if frac < 0.01 :
                convergence = True
        
    return convergence

# Plot system energy vs timestep
def Plot_SystemE(t, timestep, system_energy_vals, data):
    if t > 1000 : # PLot every 1000 timesteps
        if t % 1000 == 0 :
            # Save system energy values to datafile
            np.savetxt("SystemEnergyT100SF1.txt", system_energy_vals, delimiter = ",")
            data.to_csv('T = 100 with sf=1.csv')
            """
            # Plot total energy of system to see if system is converging
            fig, ax = plt.subplots(1,1)
            abs_energy = -1 * np.array(system_energy_vals)
            ax.plot(timestep, np.log(abs_energy))
            ax.set_xlabel('Timestep')
            ax.set_ylabel('log(System energy)')
            plt.show()
            data.to_csv('T = 100 with sf=1.csv')
            """