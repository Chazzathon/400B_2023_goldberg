# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:33:36 2023

@author: ceg30

Read MW_000.txt file to find the parameters of particles in the milky way.
"""

import numpy as np
import astropy.units as u
from pathlib import Path

def Read(filename):
    
    #opem the file
    file=open(filename,'r')
    
    #read the first line of the file to find the time. Define time units as Myr
    line1=file.readline()
    label,value=line1.split()
    time=float(value)*u.Myr
    print('Time= '+str(time))
    
    #read the second line of the file to find the total number of particles
    #described.
    line2=file.readline()
    label,value=line2.split()
    total=float(value)
    print('Number of objects= '+str(total))
    
    #close the file
    file.close()
    
    #convert all lines after line 3 to a numpy array that can be used by
    #the ParticalProperties.py script.
    data=np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    return data
    
def main():
    #define the file path to be used to open the file
    filename=Path('C:/Users/ceg30/Desktop/UA/ASTR400B/400B_2023_goldberg/Homework/Homework2/MW_000.txt')
    
    #call Read function to produce a numpy array of the data within the file
    Read(filename)
    
main()