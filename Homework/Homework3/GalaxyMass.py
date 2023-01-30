# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 13:18:40 2023

@author: ceg30
"""

import numpy as np
import astropy.units as u
from pathlib import Path
from ReadFile import Read


#This function reads from the MW_000.txt file and calculates the total mass of
#either halo, disk, or bulge mass of a galaxy.

#Inputs: filename (name of the data file to be read); particle type (integer
#number 1, 2, or 3 corresponding to the type pf particle selected)

#Returns: total_mass (total mass of halo, disk, or bulge of galaxy in units of
#10^12 solar mass)
def ComponentMass(filename,particle_type):
    
    data,time,total=Read(filename)#get data array by calling the Read function 
                                  #from ReadFile
    
    index=np.where(data['type']==particle_type)#find the indices where the
                                               #particles of chosen type are
                                               #located
    
    mass_list=list(data['m'][index])#get a numpy array of all of the mass
                                    #values using the indices found in the
                                    #previous line, the turn it into a list
    
    total_mass=sum(mass_list)*u.Msun*10**10#find the sum of the mass values
                                           #list and define the units
    
    total_mass=np.around(total_mass.to(u.Msun*10**12),3)#convert the units to
                                                        #10^12 solar masses and
                                                        #round to 3 decimal
                                                        #places
    
    return total_mass

def main():
    
    #take an input for the galaxy to be analyzed
    galaxy=input('Enter desired galaxy (MW, M33, M31):\n')
    
    #ask user for particle type to analyze
    particle_type=int(input('Enter particle type (1 for halo type, 2 for diskmatter, 3 for bulge type):\n'))
    
    #define the file path to the text file using the galaxy input
    filename=Path('C:/Users/ceg30/Desktop/UA/ASTR400B/400B_2023_goldberg/Homework/Homework3/'+galaxy+'_000.txt')
    
    #call the ComponentMass function to find the total mass of the chosen
    #particle type within  the chosen galaxy.
    total_mass=ComponentMass(filename,particle_type)
    
    print('\nMass= '+str(total_mass))
    
main()