# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 22:07:36 2022

@author: ceg30
"""
import numpy as np

def bin_xy(x,y,n_bins):
    
    max_value=max(x)
    min_value=min(x)
    value_range=max_value-min_value
    h=value_range/n_bins
    
    bins=[]
    i=min_value
    while i < max_value:
        bins.append(i)
        i+=h

    bin_array=[]
    for i2 in bins:
        bin_array.append([])
    
    
    bin_indices=np.digitize(x,bins)
    bin_array.append([])
    bins.append(i)
    
    for i in range(0,len(bin_indices)):
        bin_array[bin_indices[i]].append(y[i])
        
    for i in range(0,len(bin_array)):
        std=np.std(bin_array[i])
        mean=np.mean(bin_array[i])
        temp_list=[]
        for i2 in range(0,len(bin_array[i])):
            if abs(bin_array[i][i2]-mean)<std*3:
                    temp_list.append(bin_array[i][i2])
                
        bin_array[i]=np.median(temp_list)
        new_x=bins
        new_y=bin_array
        
    return new_x,new_y