# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:25:48 2023

@author: ceg30
"""

import numpy as np
from pathlib import Path

from ReadFile import Read

def open_file(filename,folder):
    
    path=Path(folder+'/'+filename)
    
    data=Read(path)[0]
    
    return data
    
    
folder='C:/Users/ceg30/Desktop/UA/ASTR400B/400B_2023_goldberg/ResearchProject/RA7'

mw_data=open_file('MW_801.txt',folder)
m31_data=open_file('M31_801.txt',folder)

combined=np.concatenate((mw_data,m31_data),axis=0)
length=str(len(mw_data)+len(m31_data))

header='Time      11442.9\nTotal       '+length+'\nmass in 1e10,  x,y,zin kpc and\
vx, vy, vz in km/s\ntype, m, x, y, z, vx, vy, vz'
np.savetxt('merged.txt',combined,header=header)