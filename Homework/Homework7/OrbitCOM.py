

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass import CenterOfMass

 
def vector_dif(vector1,vector2):
    ''''Function that finds the magnitude of differemce between two vectors.
    
    inputs:
        vector1: np array
            array containing the x, y, and z components of first vector for
            each snapshot
            
        vector2: np array
            array containing the x, y, and z components of second vector for
            each snapshot
            
    outputs:
        difference_mag: np array
            array containing the values for the magnitude of difference between
            the two vectors at each snapshot
        
    '''
    
    difference_mag=np.sqrt((vector1[0]-vector2[0])**2+(vector1[1]-vector2[1])**2+(vector1[2]-vector2[2])**2)
                                                      
    return difference_mag