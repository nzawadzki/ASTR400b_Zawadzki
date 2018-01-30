#Nicole Zawadzki
#ASTR400b Homework3
#Compute the mass breakdown of the Local Group

import numpy as np
import astropy.units as u
from astropy.table import Table, Column
from ReadFile import *

#create function that takes in filename and paticle type
#and outputs total mass of that type in units of 10^12 solar masses
def ComponentMass(filename, i):
    #read in data from ReadFile
    time, tot_particle, data = Read(filename)

    #index for each type of particle
    index = np.where(data['type']==i)

    #store mass data of i type as a list
    mass = data['m'][index]*1e10*u.Msun

    #sum the mass in the list in units of 10^12 Msun 
    tot_mass = np.sum(mass)/1e12

    return tot_mass

#computing halo, disk, buldge, stellar (disk+buldge), and total mass for MW
MW_halo_mass = ComponentMass('MW_000.txt', 1)
MW_disk_mass = ComponentMass('MW_000.txt', 2)
MW_bulge_mass = ComponentMass('MW_000.txt', 3)
Total_MW_stellar_mass = MW_disk_mass + MW_bulge_mass
Total_MW_mass = MW_halo_mass + MW_disk_mass + MW_bulge_mass

#computing masses for M31
M31_halo_mass = ComponentMass('M31_000.txt', 1)
M31_disk_mass = ComponentMass('M31_000.txt', 2)
M31_bulge_mass = ComponentMass('M31_000.txt', 3)
Total_M31_stellar_mass = M31_disk_mass + M31_bulge_mass
Total_M31_mass = M31_halo_mass + M31_disk_mass + M31_bulge_mass

#computing masses for M33
M33_halo_mass = ComponentMass('M33_000.txt', 1)
M33_disk_mass = ComponentMass('M33_000.txt', 2)
M33_bulge_mass = ComponentMass('M33_000.txt', 3)
Total_M33_stellar_mass = M33_disk_mass + M33_bulge_mass
Total_M33_mass = M33_halo_mass + M33_disk_mass + M33_bulge_mass

#computing total mass of local group
Total_LG_mass = Total_MW_mass + Total_M31_mass + Total_M33_mass

#compuing baryon fraction for MW, M31, M33 and local group as a whole
fbar_MW = Total_MW_stellar_mass / Total_MW_mass
fbar_M31 = Total_M31_stellar_mass / Total_M31_mass
fbar_M33 = Total_M33_stellar_mass / Total_M33_mass
fbar_LG = (Total_MW_stellar_mass + Total_M31_stellar_mass + Total_M33_stellar_mass ) / (Total_LG_mass)

#print statements for mass distribution in MW, M31, M33, and the local group
print "Total mass in MW halo (10^12Msun)", MW_halo_mass
print ""
print "Total mass in MW disk (10^12Msun)", MW_disk_mass
print ""
print "Total mass in MW bulge (10^12Msun)", MW_bulge_mass
print ""
print "Total MW mass (10^12Msun)", Total_MW_mass
print ""
print "fbar MW", fbar_MW
print ""
print "Total mass in M31 halo (10^12Msun)", M31_halo_mass
print ""
print "Total mass in M31 disk (10^12Msun)", M31_disk_mass
print ""
print "Total mass in M31 bulge (10^12Msun)", M31_bulge_mass
print ""
print "Total M31 mass (10^12Msun)", Total_M31_mass
print ""
print "fbar M31", fbar_M31
print ""
print "Total mass in M33 halo (10^12Msun)", M33_halo_mass
print ""
print "Total mass in M33 disk (10^12Msun)", M33_disk_mass
print ""
print "Total mass in M33 bulge (10^12Msun)", M33_bulge_mass
print ""
print "Total M33 mass (10^12Msun)", Total_M33_mass
print ""
print "fbar M33", fbar_M33
print ""
print "Total Local Group mass", Total_LG_mass
print ""
print "fbar local group", fbar_LG






