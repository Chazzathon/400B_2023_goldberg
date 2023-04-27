# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 09:10:23 2023

@author: ceg30

The topic of my research project is the kinematics of the merger remnant.
This code specifically looks at the question of: does the merger remnant
rotate. It will do this by aligning the merger remnant to an x,y,z coordinate
system, then plotting the velocity as a function of x.
"""

import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile
import pandas as pd
from BinningAlgorithm2 import bin_xy

''''First read in the position and velocity data of the merger remnant. Then
rotate the remnant so that it lines up with an x,y,z coordinate system'''

plt.rcParams.update({'font.size': 24})

#initialize center of mass object
COMD = CenterOfMass("M31MWFinalCombined.txt",2)

#retrieve center of mass and velocity of merger remnant
COMP = COMD.COM_P(0.1)
COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])

# Determine positions of disk particles relative to COM 
xD = COMD.x - COMP[0].value 
yD = COMD.y - COMP[1].value 
zD = COMD.z - COMP[2].value 

# total magnitude
rtot = np.sqrt(xD**2 + yD**2 + zD**2)

# Determine velocities of disk particles relatiev to COM motion
vxD = COMD.vx - COMV[0].value 
vyD = COMD.vy - COMV[1].value 
vzD = COMD.vz - COMV[2].value 

# total velocity 
vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)

# Vectors for r and v 
r = np.array([xD,yD,zD]).T # transposed 
v = np.array([vxD,vyD,vzD]).T

def RotateFrame(posI,velI):
    """a function that will rotate the position and velocity vectors
    so that the disk angular momentum is aligned with z axis. 
    
    PARAMETERS
    ----------
        posI : `array of floats`
             3D array of positions (x,y,z)
        velI : `array of floats`
             3D array of velocities (vx,vy,vz)
             
    RETURNS
    -------
        pos: `array of floats`
            rotated 3D array of positions (x,y,z) such that disk is in the XY plane
        vel: `array of floats`
            rotated 3D array of velocities (vx,vy,vz) such that disk angular momentum vector
            is in the +z direction 
    """
    
    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    # normalize the vector
    L_norm = L/np.sqrt(np.sum(L**2))


    # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
    
    # z unit vector
    z_norm = np.array([0, 0, 1])
    
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T
    
    return pos, vel

''''Function to find the velocity dispersion as a function of radius.
This is for my second question of: is the velocity dispersion curve typical
of observed ellipticals?'''

def velo_disp(V,R):
    '''
    Function to calculate the velocity dispersion of particles at increasing
    radii from the center of the mergr remnant
    
    inputs:
        V: numpy array
            array holding the magnitude of velocity for each particle in the
            remnant
            
        R: numpy array
            array holding distance of each particle from the center of the
            remnant
            
    outputs:
        r_list: list
            list containing the radii at which the velocity dispersion
            was calculated
            
        dispersion_list: list
            list containing the calculated velocity dispersion at each radius
            in r_list
            
        velo_list: list
        
            list containing the average particle velocity at each radius in the
            r_list

    '''
    
    #step size for radius
    dr=0.5
    
    #starting radius
    r=dr
    
    #maximum radius at which to analyze velocity dispersion
    R_max=100
    
    #initialize r_list, dispersion_list, and velo_list
    r_list=[]
    dispersion_list=[]
    velo_list=[]
    
    #while loop that calculates velocity dispersion of particles located
    #between r-dr and r+dr
    while r<R_max:
        
        index1=list(np.where((R<(r+dr))&(R>(r-dr))))
        V_r=V[index1]
        
        r_list.append(r)
        dispersion_list.append(np.std(V_r))
        velo_list.append(np.mean(V_r))
        
        r+=dr
    
    return np.array(r_list),np.array(dispersion_list), velo_list

def main():
    
    #get rotated galaxy components
    rn,vn=RotateFrame(r,v)
    
    #get magnitudes of positions and velocities
    V=np.sqrt(vn[:,0]**2+vn[:,1]**2+vn[:,2]**2)
    R=np.array(np.sqrt(rn[:,0]**2+rn[:,1]**2+rn[:,2]**2))
    
    '''
    rotation stuff
    '''
    
    n_bins=500
    new_r,new_v,error=bin_xy(rn[:,0],vn[:,1],n_bins)
    
    fig,ax=plt.subplots(figsize=(10,10))
    #ax.scatter(rn[:,0],vn[:,1])
    #ax.errorbar(new_r,new_v,error,fmt='.')
    ax.scatter(new_r,new_v)
    name='M31 MW Merger Remnant position and Velocity'
    ax.set(title=name, xlabel='X position (kpc)', ylabel='Y Velocity (km s$^{-1}$)')
    
    
    '''
    velocity dispersion stuff
    '''
    
    #retrieve radius and velocity dispersion of remnant
    r_list,dispersion_list,velo_list=velo_disp(V,R)
    
    #initialize matplotlib plot
    fig,ax=plt.subplots(figsize=(10,10))
    
    #plot velocity dispersion as a function of radius
    ax.scatter(r_list,dispersion_list)
    name='M31 MW Merger Remnant Velocity Dispersion'
    ax.set(title=name, xlabel='Radius (kpc)', ylabel='Velocity Dispersion (km s$^{-1}$)')
    
    #calculate v/sigma at each radius
    v_over_dv=velo_list/dispersion_list
    
    #initialize second matplotlib plot
    fig,ax2=plt.subplots(figsize=(10,10))
    
    #plot v/sigma at each radius
    ax2.scatter(r_list,v_over_dv)
    name='M31 MW Merger Remnant V/$\sigma$'
    ax2.set(title=name, xlabel='Radius (kpc)', ylabel='V/$\sigma$')
    
    #find the total average v/sigma for the entire remnant
    mean_v=np.mean(V)
    total_disp=np.std(V)
    
    
    print('Average V/sigma= '+str(np.around(mean_v/total_disp,1)))
    
main()