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
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            filename : `str`
                snapshot file
            ptype : `int; 1, 2, or 3`
                particle type to use for COM calculations
        '''
        
        ilbl= '000' + str(snap)
        ilbl=ilbl[-3:]
        self.filename='%s_'%(galaxy)+ilbl+'.txt'
        
        
        
        # read data in the given file using Read
        self.data, self.time, self.total = Read(self.filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        #self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m']
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        
        self.gname=galaxy
        
    def MassEnclosed(self, ptype,r):
        
        COM=CenterOfMass(self.filename,ptype)
        COM_p=COM.COM_P(.1)
        
        index1=np.where(self.data['type']==ptype)
        
        x_COM=COM_p[0]
        y_COM=COM_p[1]
        z_COM=COM_p[2]
        
        x_new = self.x[index1]-x_COM
        y_new = self.y[index1]-y_COM
        z_new = self.z[index1]-z_COM
        m_new=self.m[index1]
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)
        
        mass_list=[]
        
        for radius in r:
            index2 = np.where(r_new < radius)
            m2 = m_new[index2]
            mass_list.append(sum(m2))
            
        mass_array=np.array(mass_list)*u.M_sun*10**10
            
        #fig,ax=plt.subplots(figsize=(10,10))
        #ax.scatter(r,mass_array)
        #ax.set(title='Example Mass Profile for MW Halo', xlabel='Radius (kpc)', ylabel='Log(Mass Enclosed ($M_{\odot}$))')
        
        return mass_array
        
    def MassEnclosedTotal(self,r):
        
        halo_mass=self.MassEnclosed(1,r)
        disk_mass=self.MassEnclosed(2,r)
        
        if self.gname != 'M33':
            bulge_mass=self.MassEnclosed(3,r)
            total_mass=halo_mass+disk_mass+bulge_mass
            
        else:
            total_mass=halo_mass+disk_mass
            
            
        return total_mass
    
    def HernquistMass(self,radius,a,Mhalo):
        
        m_r=(Mhalo*(radius**2)/ ((a+radius)**2))#*u.Msun
        
        return m_r
    
    
    def CircularVelocity(self,ptype,r):
         
        mass_array=self.MassEnclosed(ptype,r)
        
        velocity=np.round(np.sqrt(G.to(u.kpc*u.km**2/u.s**2/u.Msun)*mass_array/r),2)*u.km/u.s
        
        return velocity
    
    def CircularVelocityTotal(self,r):
        
        total_mass_array=self.MassEnclosedTotal(r)
        
        total_velocity=np.round(np.sqrt(G.to(u.kpc*u.km**2/u.s**2/u.Msun)*total_mass_array/r),2)*u.km/u.s
        
        return total_velocity
    
    def HernquistVCirc(self,radius, a, Mhalo):
        
        VCirc=np.round(np.sqrt(G.to(u.kpc*u.km**2/u.s**2/u.Msun)*self.HernquistMass(radius,a,Mhalo)*u.Msun/radius),2)
        
        return VCirc
        

if __name__ == '__main__' : 

    r=np.arange(0.25,30.5,0.5)*u.kpc    
    
    
    MW=MassProfile('MW',0)
    MW_halop=MW.MassEnclosed(1,r)
    MW_diskp=MW.MassEnclosed(2,r)
    MW_bulgep=MW.MassEnclosed(3,r)
    MW_totalp=MW.MassEnclosedTotal(r)
    
    Mhalo_MW=max(MW_halop)/u.Msun
    
    hernquistp=[]
    a_MW=20
    for radius in r:
        MW_Hernq=MW.HernquistMass(radius/u.kpc,a_MW,Mhalo_MW)
        hernquistp.append(MW_Hernq)
    
    
    fig,ax=plt.subplots(figsize=(10,10))
    ax.semilogy(r,MW_halop,label='Halo Profile')
    ax.semilogy(r,MW_diskp,label='Disk Profile')
    ax.semilogy(r,MW_bulgep,label='Bulge Profile')
    ax.semilogy(r,MW_totalp,label='Total Profile')
    ax.semilogy(r,hernquistp,linestyle='dashed',label='Hernquist Profile')
    
    legend=ax.legend()
    ax.set(title='MW Mass Profiles', xlabel='Radius (kpc)', ylabel='Log(Mass Enclosed ($M_{\odot}$))')
    
    print('The best scale length, a, for the MW is '+str(a_MW*u.kpc))
   

    '''M31 Mass Profile'''
    
    M31=MassProfile('M31',0)
    M31_halop=M31.MassEnclosed(1,r)
    M31_diskp=M31.MassEnclosed(2,r)
    M31_bulgep=M31.MassEnclosed(3,r)
    M31_totalp=MW.MassEnclosedTotal(r)
    
    Mhalo_M31=max(M31_halop)/u.Msun
    
    hernquistp=[]
    a_M31=6
    for radius in r:
        M31_Hernq=M31.HernquistMass(radius/u.kpc,a_M31,Mhalo_M31)
        hernquistp.append(M31_Hernq)
    
    
    fig,ax=plt.subplots(figsize=(10,10))
    ax.semilogy(r,M31_halop,label='Halo Profile')
    ax.semilogy(r,M31_diskp,label='Disk Profile')
    ax.semilogy(r,M31_bulgep,label='Bulge Profile')
    ax.semilogy(r,M31_totalp,label='Total Profile')
    ax.semilogy(r,hernquistp,linestyle='dashed',label='Hernquist Profile')
    
    legend=ax.legend()
    ax.set(title='M31 Mass Profiles', xlabel='Radius (kpc)', ylabel='Log(Mass Enclosed ($M_{\odot}$))')
    
    print('The best scale length, a, for M31 is '+str(a_M31*u.kpc))
    
    
    M33=MassProfile('M33',0)
    M33_halop=M33.MassEnclosed(1,r)
    M33_diskp=M33.MassEnclosed(2,r)
    M33_totalp=MW.MassEnclosedTotal(r)
    
    Mhalo_M33=max(M33_diskp)/u.Msun
    
    hernquistp=[]
    a_M33=8
    for radius in r:
        M33_Hernq=M33.HernquistMass(radius/u.kpc,a_M33,Mhalo_M33)
        hernquistp.append(M33_Hernq)
    
    
    fig,ax=plt.subplots(figsize=(10,10))
    ax.semilogy(r,M33_halop,label='Halo Profile')
    ax.semilogy(r,M33_diskp,label='Disk Profile')
    ax.semilogy(r,M33_totalp,label='Total Profile')
    ax.semilogy(r,hernquistp,linestyle='dashed',label='Hernquist Profile',color='red')
    
    legend=ax.legend()
    ax.set(title='M33 Mass Profiles', xlabel='Radius (kpc)', ylabel='Log(Mass Enclosed ($M_{\odot}$))')
    
    print('The best scale length, a, for M33 is '+str(a_M33*u.kpc))
    
    '''Plot the rotation curves for each galaxy'''
    '''MW Rotation Curves'''
    
    MW_halov=MW.CircularVelocity(1,r)
    MW_diskv=MW.CircularVelocity(2,r)
    MW_bulgev=MW.CircularVelocity(3,r)
    MW_totv=MW.CircularVelocityTotal(r)
    
    
    hernquistv=[]
    for radius in r:
        MW_Hernq=MW.HernquistVCirc(radius,a_MW*u.kpc,Mhalo_MW)
        hernquistv.append(MW_Hernq*u.s/u.km)
    
    
    fig,ax=plt.subplots(figsize=(10,10))
    ax.scatter(r,MW_halov,label='Halo Profile')
    ax.scatter(r,MW_diskv,label='Disk Profile')
    ax.scatter(r,MW_bulgev,label='Bulge Profile')
    ax.scatter(r,MW_totv,label='Total Profile')
    ax.plot(r,hernquistv,linestyle='dashed',label='Hernquist Profile',color='red')
    
    legend=ax.legend()
    ax.set(title='MW Velocity Profiles', xlabel='Radius (kpc)', ylabel='Log(Velocity ($ms^{-1}$))')
    
    '''M31 Rotation Curves'''
    
    M31_halov=M31.CircularVelocity(1,r)
    M31_diskv=M31.CircularVelocity(2,r)
    M31_bulgev=M31.CircularVelocity(3,r)
    M31_totv=M31.CircularVelocityTotal(r)
    
    
    hernquistv=[]
    for radius in r:
        M31_Hernq=M31.HernquistVCirc(radius,a_M31*u.kpc,Mhalo_M31)
        hernquistv.append(M31_Hernq*u.s/u.km)
    
    
    fig,ax=plt.subplots(figsize=(10,10))
    ax.scatter(r,M31_halov,label='Halo Profile')
    ax.scatter(r,M31_diskv,label='Disk Profile')
    ax.scatter(r,M31_bulgev,label='Bulge Profile')
    ax.scatter(r,M31_totv,label='Total Profile')
    ax.plot(r,hernquistv,linestyle='dashed',label='Hernquist Profile',color='red')
    
    legend=ax.legend()
    ax.set(title='M31 Velocity Profiles', xlabel='Radius (kpc)', ylabel='Log(Velocity ($ms^{-1}$))')
    
    '''M33 Rotation Curves'''
    
    M33_diskv=M33.CircularVelocity(2,r)
    M33_halov=M33.CircularVelocity(1,r)
    M33_totv=M33.CircularVelocityTotal(r)
    
    
    hernquistv=[]
    for radius in r:
        M33_Hernq=M33.HernquistVCirc(radius,a_M33*u.kpc,Mhalo_M33)
        hernquistv.append(M33_Hernq*u.s/u.km)
    
    
    fig,ax=plt.subplots(figsize=(10,10))
    ax.scatter(r,M33_halov,label='Halo Profile')
    ax.scatter(r,M33_diskv,label='Disk Profile')
    ax.scatter(r,M33_totv,label='Total Profile')
    ax.plot(r,hernquistv,linestyle='dashed',label='Hernquist Profile',color='red')
    
    legend=ax.legend()
    ax.set(title='M33 Velocity Profiles', xlabel='Radius (kpc)', ylabel='Log(Velocity ($ms^{-1}$))')
    