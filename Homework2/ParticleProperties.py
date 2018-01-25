#Nicole Zawadzki
#ASTR400b Homework 2
#Using ReadFile find mass, velocity and position of particle

import numpy as np
import astropy.units as u
import math as m
from ReadFile import *

#function to calcuate 3D position, velocity, and mass
#of the (j+1)th particle of type i
def ParticleInfo(filename,i,j):
    #read in function for ReadFile
    t,nparticle,data=Read(filename)

    #index for type of particle
    index = np.where(data['type']==i)

    #defining position vectors
    x = data['x'][index]
    y = data['y'][index]
    z = data['z'][index]

    #print type(x), len(x) --test print statement

    #defining velocity vectors
    vx = data['vx'][index]
    vy = data['vy'][index]
    vz = data['vz'][index]

    #defining mass
    mass = data['m'][index]

    #calculating magnitude of distance, velocity, and mass
    distance = m.sqrt(x[j]**2 + y[j]**2 + z[j]**2)*u.kpc
    velocity = m.sqrt(vx[j]**2 + vy[j]**2 + vz[j]**2)*u.km/u.s
    Mass = mass[j]*1e10*u.Msun

    #round answers to 3 decimal points and convert distance to lightyears
    rounded_distance=np.around(distance.to(u.lyr),3)
    rounded_velocity=np.around(velocity,3)
    rounded_mass=np.around(Mass,3)

    return rounded_distance, rounded_velocity, rounded_mass

#print information for 100th particle of type 2(disk)
d,v,m = ParticleInfo('MW_000.txt',2,99)
print "3D distance=", d
print "3D velocity=", v
print "Mass=", m




