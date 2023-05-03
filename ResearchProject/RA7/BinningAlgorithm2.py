# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 13:17:30 2023

@author: ceg30

Function to bin data and find the errors for each bin.
"""

import numpy as np

def bin_xy(x,y,n_bins):
    '''
    function that takes x and y data and bins it into a set number of bins
    by averaging the y value of data points that fall within each x bin
    
    inputs:
        x: list or numpy array
            x values of plot
            
        y: list or numpy array
            y values of plot
            
        n_bins: integer
            number of desired bins
            
    returns:
        new_x: list
            binned x values to be plotted
            
        new_y: list
            binned y values to be plotted
            
        error: list
            error for each binned y value calculated from standard deviation
            of all original y values in each given bin
    '''
    
    
    #find max and minimum bin values as well as spacing in between each bin
    max_value=max(x)
    min_value=min(x)
    value_range=max_value-min_value
    h=value_range/n_bins
    
    #create list of the bins
    bins=[]
    i=min_value
    while i < max_value:
        bins.append(i)
        i+=h
        
    #create an array of empty lists that will hold the y values that fall into
    #each bin
    bin_array=[]
    for i2 in bins:
        bin_array.append([])
    
    #create a list of indices that correspond to the bin that each xy pair
    #belongs to. Each y value will be placed into its corresponding assigned
    #index/list in the bin_array
    bin_indices=np.digitize(x,bins)
    
    bin_array.append([])
    bins.append(i)
    
    #add each y value into the correct bin list within the bin_array
    for i in range(0,len(bin_indices)):
        bin_array[bin_indices[i]].append(y[i])
    
    #loop through the lists within the bin_arrayp, finding the standard
    #deviation and mean of the y values within each bin
    error=[]
    for i in range(0,len(bin_array)):
        std=np.std(bin_array[i])
        mean=np.mean(bin_array[i])
        
        #loop through the individual y values within each bin and check that
        #they are within 3 standard deviations of the mean y value
        temp_list=[]
        for i2 in range(0,len(bin_array[i])):
            
            #Add the accepted values to a new list, then take the mean as the y
            #value for that bin
            if abs(bin_array[i][i2]-mean)<std*3:
                    temp_list.append(bin_array[i][i2])
         
        #replace the list of y values in the current bin with the final y value
        #and set the error as the standard deviation of the sigma clipped list
        bin_array[i]=np.median(temp_list)
        error.append(np.std(temp_list))
        
    new_x=bins
    new_y=bin_array
        
    return new_x,new_y,error