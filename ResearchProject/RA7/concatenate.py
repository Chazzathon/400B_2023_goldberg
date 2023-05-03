# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:25:48 2023

@author: ceg30

File to sitch together final M31 and MW files so that all of their particles
can be analyzed as one galaxy.
"""

#import modules
import numpy as np
from pathlib import Path
from ReadFile import Read

def open_file(filename,folder):
    '''
    function that reads the data from a file into a np array

    Parameters
    ----------
    filename : str
        string of the name of the file to be opened
        
    folder : str
        string of the file path to the folder where the file is located

    Returns
    -------
    data : nparray
        numpy array with all of the data from the file that was read

    '''
    
    #create a path to the file using the folder and filename
    path=Path(folder+'/'+filename)
    
    #create the numpy array by reading in the file using the Read function from
    #ReadFile
    data=Read(path)[0]
    
    return data
    
#specify the folder that both files are in  
folder='C:/Users/ceg30/Desktop/UA/ASTR400B/400B_2023_goldberg/ResearchProject/RA7'

#create np arrays of both files using the open_file function
mw_data=open_file('MW_801.txt',folder)
m31_data=open_file('M31_801.txt',folder)

#combine the data from both files into one object
combined=np.concatenate((mw_data,m31_data),axis=0)
length=str(len(mw_data)+len(m31_data))

#specify the header of the file that is about to be created
header='Time      11442.9\nTotal       '+length+'\nmass in 1e10,  x,y,zin kpc and\
vx, vy, vz in km/s\ntype, m, x, y, z, vx, vy, vz'

#save the merged data to a new txt file
np.savetxt('merged.txt',combined,header=header)