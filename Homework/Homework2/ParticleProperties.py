# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:47:58 2023

@author: ceg30
"""

import numpy as np
import astropy.units as u
import ReadFile
from pathlib import Path

#This function reads from the MW_000.txt file and calculates the distance,
#velocity, and mass of the chosen particle.

#Inputs: filename (name of the data file to be read); particle type (integer
#number 1, 2, or 3 corresponding to the type pf particle selected);
#particle_num (the index corresponding to the ith particle of type 1,2 or 3.)

#Returns: distance in light years, velocity in km/s, and mass in Msun
def ParticleInfo(filename,particle_type,particle_num):
    
    data=ReadFile.Read(filename)#get data array by calling the Read function 
                                #from ReadFile
    
    index=np.where(data['type']==particle_type)#find the indices where the
                                               #particles of chosen type are
                                               #located
    
    #filter out the other particle types using the indices just found. Then
    #select the x,y, and z locations of the chosen particle, identified by the
    #"particle_num" index.
    x=data['x'][index][particle_num]
    y=data['y'][index][particle_num]
    z=data['z'][index][particle_num]
    
    #filter out the other particle types using the indices just found. Then
    #select the x,y, and z velocities of the chosen particle, identified by the
    #"particle_num" index.
    vx=data['vx'][index][particle_num]
    vy=data['vy'][index][particle_num]
    vz=data['vz'][index][particle_num]
    
    #filter out the other particle types using the indices just found. Then
    #select the mass value of the chosen particle, identified by the
    #"particle_num" index.
    mass=data['m'][index][particle_num]
    
    #find the magnitude of the distance by taking the square root of the sum of
    #squares of the x, y, and z locations. Then define the units as kpc and
    #convert to light years using astropy units. Round to 3 decimal places.
    distance=(np.sqrt(x**2+y**2+z**2))*u.kpc
    distance=np.around(distance.to(u.lyr),3)
    
    #find the magnitude of the velocity by taking the square root of the sum of
    #squares of the x, y, and z locations. Then define the units as km/s using
    #astropy units. Round to 3 decimal places.
    velocity=(np.sqrt(vx**2+vy**2+vz**2))*u.km/u.s
    velocity=np.around(velocity,3)
    
    #define the mass units as 1E10 solar masses, then convert to solar masses.
    #Round to 3 decimal places.
    mass=mass*(u.M_sun*10**10)
    mass=np.around(mass.to(u.M_sun),3)
    
    return distance, velocity, mass

def main():
    
    #define the file path to be used to open the file
    filename=Path('C:/Users/ceg30/Desktop/UA/ASTR400B/400B_2023_goldberg/Homework/Homework2/MW_000.txt')
    
    #ask user for particle type to search for
    particle_type=int(input('Enter particle type (1 for dark matter, 2 for disk stars, 3 for bulge stars):\n'))
    
    #ask user for desired particle number of the chosen type. Subtract by 1 to
    #account for python's zero based indexing(i.e. particle #1 becomes index 0)
    particle_num=int(input('Enter index number of desired particle type:\n'))
    particle_num-=1
    
    #call ParticleInfo function and retrieve the parameters it returns
    distance, velocity, mass=ParticleInfo(filename,particle_type,particle_num)
    
    #print returned values from ParticleInfo function
    print('\nDistance= '+str(distance))
    print('Velocity= '+str(velocity))
    print('Mass= '+str(mass))
    
main()