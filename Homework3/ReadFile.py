#Nicole Zawadzki
#ASTR400B Homework 2
#Read in file and return information

import numpy as np
import astropy.units as u

filename='MW_000.txt'

#define function that takes name of file as input
def Read(filename):
    f = open(filename, 'r')
    #read in first line and store time in units of 10Myr
    line1 = f.readline()
    label1,value1=line1.split()
    time=float(value1)*1e10*u.Myr
    #print 'Time:', time
    
    #read the second line and store total number of particles
    line2=f.readline()
    label2,value2=line2.split()
    tot_particle=float(value2)
    #print 'Total Particles=', tot_particle

    f.close()
    
    #store rest of file as a data array
    data=np.genfromtxt(filename, dtype=None, names=True,skip_header=3)

    return time, tot_particle, data

t,nparticle,d=Read('MW_000.txt')

                                           
