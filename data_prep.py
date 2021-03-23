# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 10:54:25 2020

@author: s1708916
"""
# Should you make this a class?
import math
import numpy as np
import pandas as pd

def DegreeToDistance(data):
    r_earth = 6371.0 #km
    data['x'] = data['Longitude'] * math.pi * r_earth / 180
    data['y'] = data['Latitude'] * math.pi * r_earth / 180
    return data

def CentredDistance(data):
    min_x = data['x'].min()
    min_y = data['y'].min()
    data['x'] = data['x'] - min_x
    data['y'] = data['y'] - min_y
    data.sort_values(by = 'x', inplace = True, ascending = True)
    return data
    

# Coarse data
def CoarseData(data, resolution):
    res = resolution
    
    x_max = data['x'].max()
    x_min = data['x'].min()
    y_max = data['y'].max()
    y_min = data['y'].min()
    
    i_range = math.ceil((x_max - x_min)/res)
    j_range = math.ceil((y_max - y_min)/res)
    
    coarse_map = pd.DataFrame(columns = ['Longitude', 'Latitude', 'Population', 'x', 'y'])
    
    x_vals = []
    y_vals = []
    pop_vals = []
    long = []
    lat = []
    
    for i in range(i_range): #Loop over x vals
        for j in range(j_range): # Loop over y vals
            grid = data[(data['x'] > x_min + res*i) & (data['x'] < (x_min + res*(i+1)))]
            grid = grid[(grid['y'] > y_min + res*j) & (grid['y'] < (y_min + res*(j+1)))]
            
            x_vals.append(grid['x'].mean())
            y_vals.append(grid['y'].mean())
            pop_vals.append(grid['Population'].mean()) # should this not be sum
            long.append(grid['Longitude'].mean())
            lat.append(grid['Latitude'].mean())
            
            
    coarse_map['x'] = x_vals
    coarse_map['y'] = y_vals
    coarse_map['Population'] = pop_vals
    coarse_map['Longitude'] = long
    coarse_map['Latitude'] = lat
   
    # Remove zero population values
    coarse_map = coarse_map.dropna()
    coarse_map = coarse_map.reset_index(drop = True)
    
    # Should technically change to population density but wouldn't change that much
    return coarse_map
            
#Reduced data
def ReducedData(data):
    
    test_point = data.sample(1)
    
    test_point = test_point.index
    
    print(test_point)
    
    x = data['x'].loc[test_point]
    y = data['y'].loc[test_point]
    
    red_grid = data[(data['x'] < x)]
    
    #& (data['x'] > (x - 50))]
    #red_grid = red_grid[(red_grid['y'] < (y + 50)) & (red_grid['y'] > (y - 50))]
    
    return red_grid
    
    