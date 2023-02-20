

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


def OrbitCOM(galaxy,start,end,n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
        galaxy: str
            string with the name of the galaxy to be analyzed
            
        start: int
            number of first snapshot file to analyze
            
        end: int
            number of final snapshot file to analyze
            
        n: int
            interval between files that are analyzed
          
    outputs:
        .txt file
            File containing the center of mass position and velocity vectors
            for the chosen galaxy at each snapshot analyzed
    """
    
    # compose the filename for output
    fileout='Orbit'+str(galaxy)+'.txt'
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta=0.1
    if galaxy=='M33':
        volDEC=4
    else:
        volDEC=2

    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids=np.arange(start,end,n)
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit=np.zeros([len(snap_ids),7])
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        ilbl='000'+str(snap_id)
        ilbl=ilbl[-3:]
        filename='%s'%(galaxy)+'_'+ilbl+'.txt'
        
        # Initialize an instance of CenterOfMass class, using disk particles
        COM=CenterOfMass(filename,ptype=2)
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COM_pos=COM.COM_P(delta,volDEC)
        x_com=COM_pos[0]
        y_com=COM_pos[1]
        z_com=COM_pos[2]
        
        COM_vel=COM.COM_V(x_com,y_com,z_com)
        
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        orbit[i][0]=COM.time.value/1000  #add time to 1st column
        orbit[i][1]=COM_pos[0].value#add x pos to second column
        orbit[i][2]=COM_pos[1].value#add y pos to third column
        orbit[i][3]=COM_pos[2].value#add z pos to 4th column
        orbit[i][4]=COM_vel[0].value#add x vel to 5th column
        orbit[i][5]=COM_vel[1].value#add y vel to 6th column
        orbit[i][6]=COM_vel[2].value#add z vel to 7th column
        
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        



# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 

''''Create .txt file with the center of mass position and velocity vectors
for the chosen galaxy at each snapshot analyzed. Commented out once file has
been created for each galaxy.'''
#OrbitCOM('M33',0,800,5)


# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
MW_file=file=open('OrbitMW.txt','r')
MW_data=np.genfromtxt(MW_file,dtype=None,names=True,skip_header=0)


M31_file=file=open('OrbitM31.txt','r')
M31_data=np.genfromtxt(M31_file,dtype=None,names=True,skip_header=0)

M33_file=file=open('OrbitM33.txt','r')
M33_data=np.genfromtxt(M33_file,dtype=None,names=True,skip_header=0)


# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
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



# Determine the magnitude of the relative position and velocities 

# of MW and M31
MW_M31_dr=vector_dif([MW_data['x'],MW_data['y'],MW_data['z']],[M31_data['x'],M31_data['y'],M31_data['z']])
MW_M31_dv=vector_dif([MW_data['vx'],MW_data['vy'],MW_data['vz']],[M31_data['vx'],M31_data['vy'],M31_data['vz']])

# of M33 and M31
M33_M31_dr=vector_dif([M33_data['x'],M33_data['y'],M33_data['z']],[M31_data['x'],M31_data['y'],M31_data['z']])
M33_M31_dv=vector_dif([M33_data['vx'],M33_data['vy'],M33_data['vz']],[M31_data['vx'],M31_data['vy'],M31_data['vz']])


# Plot the Orbit of the galaxies 
#################################
fig,ax=plt.subplots(1,2,figsize=(10,5))

ax[0].set(title='MW and M31 Relative Distance', xlabel='Time (Gyr)', ylabel='Distance (kpc)')
ax[0].plot(MW_data['t'],MW_M31_dr)
print(len(MW_data['t']))

ax[1].set(title='M33 and M31 Relative Distance', xlabel='Time (Gyr)', ylabel='Distance (kpc)')
ax[1].plot(MW_data['t'],M33_M31_dr)



# Plot the orbital velocities of the galaxies 
#################################
fig,ax2=plt.subplots(1,2,figsize=(10,5))

ax2[0].set(title='MW and M31 Relative Velocity', xlabel='Time (Gyr)', ylabel='Velocity (km s$^{-1}$)')
ax2[0].plot(MW_data['t'],MW_M31_dv)

ax2[1].set(title='M33 and M31 Relative Velocity', xlabel='Time (Gyr)', ylabel='Velocity (km s$^{-1}$)')
ax2[1].plot(MW_data['t'],M33_M31_dv)

''''Answer Questions'''

#1
print('1: The MW and M31 will have two close encounters before merging.\n')

#2
print('2: When the galaxies are far apart their relative velocities are small, \
but they speed up as they approach each other.\n')
      
#3
print('3: The MW and M31 merge at about 6 Gyrs. The radius of M33s orbit becomes \
smaller after the merger.\n')
    
#bonus:
print('4: M33s orbit is decaying at ~15 Kpc/Gyr. So once it is at 75 Kpc, \
it will take 5 Gyr for M33 to merge with the rest of the remnant.\n')
