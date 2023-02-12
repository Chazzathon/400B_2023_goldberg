# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 19:45:23 2023

@author: ceg30
"""

import numpy as np
import astropy.units as u
import astropy.table as tbl
from astropy.constants import G

from ReadFile import Read
import matplotlib.pyplot as plt

from CenterOfMass import CenterOfMass

class MassProfile:
    
    def __init__(self, galaxy, snap):
        ''' Class to calculate the mass distribution and velocity profiles of
            galaxies
            
            PARAMETERS
            ----------
            galaxy : `str`
                galaxy to be analyzed
            snap : `int`
                specified time in simulation at which galaxy is analyzed
        '''
        
        #add a string of the file number to the value "000"
        ilbl= '000' + str(snap)
        #remove all but the last 3 digits
        ilbl=ilbl[-3:]
        self.filename='%s_'%(galaxy)+ilbl+'.txt'
        
        
        
        # read data in the given file using Read
        self.data, self.time, self.total = Read(self.filename)                                                                                             

        #create an array to store indexes of particles                                
        #self.index = np.where(self.data['type'] == ptype)

        # store the mass and positions of the particles
        self.m = self.data['m']
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        
        #define the name of the galaxy for future use
        self.gname=galaxy
        
    def MassEnclosed(self, ptype,r):
        ''''function to find the mass enclosed within different radii
        
        inputs:
            ptype: 'int'
            specifies the type of particle to be used to calculate the mass from
            
            r: 'np array'
            array holding the r values at which the enclosed mass is to be
            calculated
            
        outputs:
            mass_array: 'np array'
            the masses enclosed wihthin the radii specified by r.
        
        '''
        
        #find the center of mass of the desired particle type
        COM=CenterOfMass(self.filename,ptype)
        COM_p=COM.COM_P(.1)
        
        #specify the indices of data at which the desired particle type is stored
        index1=np.where(self.data['type']==ptype)
        
        #define the x, y, and z center of masses
        x_COM=COM_p[0]
        y_COM=COM_p[1]
        z_COM=COM_p[2]
        
        #find the positions of the particles relative to the center of mass.
        #THis is then used to find the radius of each particle from the center
        #of the galaxy
        x_new = self.x[index1]-x_COM
        y_new = self.y[index1]-y_COM
        z_new = self.z[index1]-z_COM
        m_new=self.m[index1]
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)
        
        #define list that will hold the enclosd mass values
        mass_list=[]
        
        #loop through the radii in the r list. Select particles that are
        #located within the given radius, and add their combined mass to the
        #mass list
        for radius in r:
            index2 = np.where(r_new < radius)
            m2 = m_new[index2]
            mass_list.append(sum(m2))
         
        #turn the mass list into a numpy array
        mass_array=np.array(mass_list)*u.M_sun*10**10
            
        
        return mass_array
        
    def MassEnclosedTotal(self,r):
        ''''
        functuion to find the total mass enclosed by all particle types at
        given radii
        
        inputs:
            r: 'np array'
            array holding the r values at which the enclosed mass is to be
            calculated
            
        outputs: 
            total_mass: 'np array'
            the masses enclosed wihthin the radii specified by r.
        '''
        
        #find the arrays of halo and disk mass enclosed by radii defined by r
        halo_mass=self.MassEnclosed(1,r)
        disk_mass=self.MassEnclosed(2,r)
        
        #m33 does not have a bulge, so if the specified galaxy is not M33 the
        #mass array for the bulge is found. All 3 arrays are then added
        #together to find the total mass enclosed by radii defined by r.
        if self.gname != 'M33':
            bulge_mass=self.MassEnclosed(3,r)
            total_mass=halo_mass+disk_mass+bulge_mass
          
        #if the galaxy is M33, the total mass is just the halo plus disk mass
        else:
            total_mass=halo_mass+disk_mass
            
            
        return total_mass
    
    def HernquistMass(self,radius,a,Mhalo):
        '''
        function to find the mass enclosed within a radius based upon the 
        hernquist profile
        
        inputs:
            radius: 'astropy quantity' in units of kpc
            distance from center of galaxy where m enclosed is to be evaluated
            
            a: 'astropy quantity' in units of kpc
            scale factor of galaxy
            
        outputs:
            m_r: 'astropy quantity' in units of Msun
            mass enclosed by given radius
        

        '''
        
        #calculate the mass enclosed by the radius using the hernquist profile
        #m=M_halo * r^2 / (a+r)^2
        m_r=(Mhalo*(radius**2)/ ((a+radius)**2))
        
        return m_r
    
    
    def CircularVelocity(self,ptype,r):
        '''
        function to calculate the orbital velocity of objects at their radius
        
        inputs:
            ptype: 'int'
            specifies the type of particle to be used to calculate the velocity for
            
            r: 'np array'
            array holding the r values at which the enclosed mass is to be
            calculated
            
        outputs:
            velocity: numpy array in units of km/s
            array holding the orbital velocities of objects
        
        '''
         
        #get the enclosed mass at specified radii due to objects of the desired
        #type
        mass_array=self.MassEnclosed(ptype,r)
        
        #calculate the velocity of each object using the orbital velocity equation
        #v=sqrt(GM/r)
        velocity=np.round(np.sqrt(G.to(u.kpc*u.km**2/u.s**2/u.Msun)*mass_array/r),2)*u.km/u.s
        
        return velocity
    
    def CircularVelocityTotal(self,r):
        '''
        function to calculate the orbital velocity of objects at their radius
        using the total enclosed mass due to all particle types
        
        inputs:
            r: 'np array'
            array holding the r values at which the enclosed mass is to be
            calculated
            
        outputs:
            total_velocity: numpy array in units of km/s
            array holding the orbital velocities of objects
        
        '''
        
        #retrieve the total mass enclosed per radius
        total_mass_array=self.MassEnclosedTotal(r)
        
        #calculate the velocity of each object using the orbital velocity equation
        #v=sqrt(GM/r)
        total_velocity=np.round(np.sqrt(G.to(u.kpc*u.km**2/u.s**2/u.Msun)*total_mass_array/r),2)*u.km/u.s
        
        return total_velocity
    
    def HernquistVCirc(self,radius, a, Mhalo):
        '''
        function to calculate the orbital velocity of objects using an enclosed
        mass specified by the hernquist profile
        '''

        #calculate the velocity of each object using the orbital velocity equation
        #v=sqrt(GM/r)
        VCirc=np.round(np.sqrt(G.to(u.kpc*u.km**2/u.s**2/u.Msun)*self.HernquistMass(radius,a,Mhalo)*u.Msun/(radius)),2)
        
        return VCirc
        

if __name__ == '__main__' : 

    #define values of r at which to eveluate enclosed mass and velocities
    r=np.arange(0.25,30.5,0.5)*u.kpc    
    
    '''MW Mass Profile'''
    
    #define center of mass object for milky way
    MW=MassProfile('MW',0)
    
    #find the mass profiles of MW halo, disk and bulge
    MW_halop=MW.MassEnclosed(1,r)
    MW_diskp=MW.MassEnclosed(2,r)
    MW_bulgep=MW.MassEnclosed(3,r)
    MW_totalp=MW.MassEnclosedTotal(r)
    
    #find the mass of the halo
    Mhalo_MW=max(MW_halop)/u.Msun
    
    #create empty list to store enclosed mass values derived from the hernquist
    #profile
    hernquistp=[]
    
    #specify scale factor of MW
    a_MW=30
    
    #find the enclosed mass according to hernquist profile at each radius
    for radius in r:
        MW_Hernq=MW.HernquistMass(radius/u.kpc,a_MW,Mhalo_MW)
        hernquistp.append(MW_Hernq)
    
    #set up graph and plot MW mass profiles on log graphs
    fig,ax=plt.subplots(figsize=(10,10))
    ax.semilogy(r,MW_halop,label='Halo Profile')
    ax.semilogy(r,MW_diskp,label='Disk Profile')
    ax.semilogy(r,MW_bulgep,label='Bulge Profile')
    ax.semilogy(r,MW_totalp,label='Total Profile')
    ax.semilogy(r,hernquistp,linestyle='dashed',label='Hernquist Profile')
    
    #create legend and label graph axes
    legend=ax.legend()
    ax.set(title='MW Mass Profiles', xlabel='Radius (kpc)', ylabel='Log(Mass Enclosed ($M_{\odot}$))')
    
    print('The best scale length, a, for the MW is '+str(a_MW*u.kpc))
   

    '''M31 Mass Profile'''
    
    #define center of mass object for m31
    M31=MassProfile('M31',0)
    
    #find the mass profiles of M31 halo, disk and bulge
    M31_halop=M31.MassEnclosed(1,r)
    M31_diskp=M31.MassEnclosed(2,r)
    M31_bulgep=M31.MassEnclosed(3,r)
    M31_totalp=MW.MassEnclosedTotal(r)
    
    #find the mass of the halo
    Mhalo_M31=max(M31_halop)/u.Msun
    
    #create empty list to store enclosed mass values derived from the hernquist
    #profile
    hernquistp=[]
    
    #specify scale factor of M31
    a_M31=15
    
    #find the enclosed mass according to hernquist profile at each radius
    for radius in r:
        M31_Hernq=M31.HernquistMass(radius/u.kpc,a_M31,Mhalo_M31)
        hernquistp.append(M31_Hernq)
    
    #set up graph and plot M31 mass profiles on log graphs
    fig,ax=plt.subplots(figsize=(10,10))
    ax.semilogy(r,M31_halop,label='Halo Profile')
    ax.semilogy(r,M31_diskp,label='Disk Profile')
    ax.semilogy(r,M31_bulgep,label='Bulge Profile')
    ax.semilogy(r,M31_totalp,label='Total Profile')
    ax.semilogy(r,hernquistp,linestyle='dashed',label='Hernquist Profile')
    
    #create legend and label graph axes
    legend=ax.legend()
    ax.set(title='M31 Mass Profiles', xlabel='Radius (kpc)', ylabel='Log(Mass Enclosed ($M_{\odot}$))')
    
    print('The best scale length, a, for M31 is '+str(a_M31*u.kpc))
    
    '''M33 Mass Profile'''
    
    #define center of mass object for m31
    M33=MassProfile('M33',0)
    
    #find the mass profiles of M31 halo and disk
    M33_halop=M33.MassEnclosed(1,r)
    M33_diskp=M33.MassEnclosed(2,r)
    M33_totalp=MW.MassEnclosedTotal(r)
    
    #find the mass of the halo
    Mhalo_M33=max(M33_halop)/u.Msun
    
    #create empty list to store enclosed mass values derived from the hernquist
    #profile
    hernquistp=[]
    
    #specify scale factor of M33
    a_M33=8
    
    #find the enclosed mass according to hernquist profile at each radius
    for radius in r:
        M33_Hernq=M33.HernquistMass(radius/u.kpc,a_M33,Mhalo_M33)
        hernquistp.append(M33_Hernq)
    
    #set up graph and plot M33 mass profiles on log graphs
    fig,ax=plt.subplots(figsize=(10,10))
    ax.semilogy(r,M33_halop,label='Halo Profile')
    ax.semilogy(r,M33_diskp,label='Disk Profile')
    ax.semilogy(r,M33_totalp,label='Total Profile')
    ax.semilogy(r,hernquistp,linestyle='dashed',label='Hernquist Profile',color='red')
    
    #create legend and label graph axes
    legend=ax.legend()
    ax.set(title='M33 Mass Profiles', xlabel='Radius (kpc)', ylabel='Log(Mass Enclosed ($M_{\odot}$))')
    
    print('The best scale length, a, for M33 is '+str(a_M33*u.kpc))
    
    '''Plot the rotation curves for each galaxy'''
    '''MW Rotation Curves'''
    
    #find the velocity profiles of MW halo, disk and bulge
    MW_halov=MW.CircularVelocity(1,r)
    MW_diskv=MW.CircularVelocity(2,r)
    MW_bulgev=MW.CircularVelocity(3,r)
    MW_totv=MW.CircularVelocityTotal(r)
    
    #create empty list to store velocity values derived from the hernquist
    #profile
    hernquistv=[]
    
    #find the enclosed mass according to hernquist profile at each radius
    for radius in r:
        MW_Hernq=MW.HernquistVCirc(radius,a_MW*u.kpc,Mhalo_MW)
        hernquistv.append(MW_Hernq*u.s/u.km)
    
    #set up graph and plot MW velocity profiles on log graphs
    fig,ax=plt.subplots(figsize=(10,10))
    ax.scatter(r,MW_halov,label='Halo Profile')
    ax.scatter(r,MW_diskv,label='Disk Profile')
    ax.scatter(r,MW_bulgev,label='Bulge Profile')
    ax.scatter(r,MW_totv,label='Total Profile')
    ax.plot(r,hernquistv,linestyle='dashed',label='Hernquist Profile',color='red')
    
    #create legend and label graph axes
    legend=ax.legend()
    ax.set(title='MW Velocity Profiles', xlabel='Radius (kpc)', ylabel='Log(Velocity ($ms^{-1}$))')
    
    '''M31 Rotation Curves'''
    
    #find the velocity profiles of M31 halo, disk and bulge
    M31_halov=M31.CircularVelocity(1,r)
    M31_diskv=M31.CircularVelocity(2,r)
    M31_bulgev=M31.CircularVelocity(3,r)
    M31_totv=M31.CircularVelocityTotal(r)
    
    #create empty list to store velocity values derived from the hernquist
    #profile
    hernquistv=[]
    
    #find the enclosed mass according to hernquist profile at each radius
    for radius in r:
        M31_Hernq=M31.HernquistVCirc(radius,a_M31*u.kpc,Mhalo_M31)
        hernquistv.append(M31_Hernq*u.s/u.km)
    
    #set up graph and plot M31 velocity profiles on log graphs
    fig,ax=plt.subplots(figsize=(10,10))
    ax.scatter(r,M31_halov,label='Halo Profile')
    ax.scatter(r,M31_diskv,label='Disk Profile')
    ax.scatter(r,M31_bulgev,label='Bulge Profile')
    ax.scatter(r,M31_totv,label='Total Profile')
    ax.plot(r,hernquistv,linestyle='dashed',label='Hernquist Profile',color='red')
    
    #create legend and label graph axes
    legend=ax.legend()
    ax.set(title='M31 Velocity Profiles', xlabel='Radius (kpc)', ylabel='Log(Velocity ($ms^{-1}$))')
    
    '''M33 Rotation Curves'''
    
    #find the velocity profiles of M33 halo and disk
    M33_diskv=M33.CircularVelocity(2,r)
    M33_halov=M33.CircularVelocity(1,r)
    M33_totv=M33.CircularVelocityTotal(r)
    
    #create empty list to store velocity values derived from the hernquist
    #profile
    hernquistv=[]
    
    #find the enclosed mass according to hernquist profile at each radius
    for radius in r:
        M33_Hernq=M33.HernquistVCirc(radius,a_M33*u.kpc,Mhalo_M33)
        hernquistv.append(M33_Hernq*u.s/u.km)
    
    #set up graph and plot M31 velocity profiles on log graphs
    fig,ax=plt.subplots(figsize=(10,10))
    ax.scatter(r,M33_halov,label='Halo Profile')
    ax.scatter(r,M33_diskv,label='Disk Profile')
    ax.scatter(r,M33_totv,label='Total Profile')
    ax.plot(r,hernquistv,linestyle='dashed',label='Hernquist Profile',color='red')
    
    #create legend and label graph axes
    legend=ax.legend()
    ax.set(title='M33 Velocity Profiles', xlabel='Radius (kpc)', ylabel='Log(Velocity ($ms^{-1}$))')
    