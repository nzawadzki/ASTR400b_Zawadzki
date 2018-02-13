#Nicole Zawadzki
#ASTR 400B Homework 5
#Determine the mass distribution of each galaxy at a certain snap shot and use this to determine each galaxy's rotation curve
#For plotting: All plots are saved in this folder.
#To recreate a plot with this code, uncomment the section of interest
import numpy as np
import astropy.units as u
#import necessary files
from ReadFile import Read
from CenterOfMass import CenterOfMass

#import relavent plotting modules
import matplotlib.pyplot as plt
#get_ipython().magic('matplotlib inline')

#constants
import astropy.constants as constant

#create a class
class MassProfile:
    
    def __init__(self, galaxy, snap):
        #use inputs to reconstruct filename
        #add string of the filenumber to find value "000"
        ilbl = '000' + str(snap)
        #remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl + '.txt'

        #read in data from ReadFile
        self.time, self.total, self.data = Read(self.filename)
        
        #store position and mass data from ReadFile
        self.x = self.data['x']
        self.y = self.data['y']
        self.z = self.data['z']
        self.m = self.data['m']
        
        #store name of the galaxy as a global property
        self.gname = galaxy
    
    #function to compute the mass enclosed within a given radii of the center  of mass
    #takes particle type and array of radii as input and returns array of masses in solar masses
    def MassEnclosed(self, ptype, r):
        marray = []
        #find the center of mass position using function from CenterOfMass that takes in filename and particle type and returns x,y,z components of center of mass to a certain toleance
        COM = CenterOfMass(self.filename, ptype)
        COMP = COM.COM_P(1.0)
        
        #array if indixes of partices of desired type
        index = np.where(self.data['type'] == ptype)
        #find postiion of all particles relative to center of mass
        m2 = self.m[index]
        x2 = self.x[index] - COMP[0]
        y2 = self.y[index] - COMP[1]
        z2 = self.z[index] - COMP[2]
        r2 = np.sqrt(x2**2 + y2**2 + z2**2)
        #loop over radius array to define particles that are enclosed within the radius given at each element
        for radius in r:
            #find and store particle masses that are within given radius
            index2 = np.where(r2 < radius)
            mnew = m2[index2]
            #sum to find mass enclosed
            msum = np.sum(mnew)
            #append to mass array
            marray.append(msum)
            
        return marray*u.Msun
    
    #takes in array of radii and return array of total masses at each radius
    def MassEnclosedTotal(self,r):
        #initialize total mass array
        mtotarray = []

        #use MassEnclosed to find the mass for each partice type
        HaloMassEnclosed = self.MassEnclosed(1,r)
        DiskMassEnclosed = self.MassEnclosed(2,r)
        #caveat to account for M33 having no buldge
        if self.gname != 'M33':
            BuldgeMassEnclosed = self.MassEnclosed(3,r)
            #looping over radius element
            for element in range(len(r)):
                mtotal = (HaloMassEnclosed[element] + DiskMassEnclosed[element] + BuldgeMassEnclosed[element])
                mtotarray.append(mtotal)
        else:
            #do the same thing but only sum halo and disk
            for element in range(len(r)):
                mtotal = (HaloMassEnclosed[element] + DiskMassEnclosed[element])
                mtotarray.append(mtotal)

        return mtotarray*u.Msun 
    
    #function that takes radius, a scale factor, and the halo mass to compute the mass enlosed awithin a radius using a theoretical profile
    def HernquistMass(self, r, a, Mhalo):
        #mass profile
        M_r = Mhalo*r**2 / (a+r)**2
        
        return M_r*u.Msun
    
    #calculate the circular speed at each radius using the mass enclosed
    #takes particle type and radius array as input and returns circular speed rounded to two decimal places
    def CircularVelocity(self, ptype, r):
        #initialize circular velocity array
        V = []
        #force of gravity equals centripital accelaration
        #solve for velocity
        G = constant.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        #Fg = G*Mm/r**2=mv**2/r
        #loop over radius array
        #define counter
        count = 0
        for radius in r:
            v = np.sqrt(G*self.MassEnclosed(ptype,r)[count]/(radius*u.kpc))
            V.append(v)
            count = count + 1
        #print V
        return V
    
    #function to calculate the total circular velocity of (buldge+disk+halo) at each input radii
    def CircularVelocityTotal(self,r):
        #initializing total circular velocity array
        Vtot = []
        #circular velocity of each component
        vHalo = self.CircularVelocity(1,r)
        vDisk = self.CircularVelocity(2,r)
        #caveat for M33
        if self.gname != 'M33':
            vBuldge = self.CircularVelocity(3,r)
            for element in range(len(r)):
                Vtot.append(vHalo[element] + vDisk[element] + vBuldge[element])
        else:
            for element in range(len(r)):
                Vtot.append(vHalo[element] + vDisk[element])
               
        return Vtot

    #function that computes circular speed using the Hernquist mass profile
    #takes radius, scale factor (a), and halo mass as input, and returns circular speed to two decimal places
    def HernquistVCirc(self,r, a, Mhalo):
        HernquistV = []
        G = constant.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        M = self.HernquistMass(r, a, Mhalo)*u.Msun
        count = 0
        for radius in r:
            vhernquist = np.sqrt(G*M[count]/(radius*u.kpc))
            HernquistV.append(vhernquist)
            count = count+1
        
        return HernquistV
 
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
#radius array out to 30kpc
rarr = np.logspace(-1,1.477,10)
print len(rarr)
#Create class objects for each galaxy
MWProfile = MassProfile('MW',0)
M31Profile = MassProfile('M31',0)
M33Profile = MassProfile('M33',0)

#Mass Components of MW
MWHaloMass = MWProfile.MassEnclosed(1,rarr)
print type(MWHaloMass)
MWDiskMass = MWProfile.MassEnclosed(2,rarr)
MWBuldgeMass = MWProfile.MassEnclosed(3,rarr)
MWTotMass = MWProfile.MassEnclosedTotal(rarr)
#print len(MWHaloMass)
#Hernquist Mass Profile
print MWHaloMass[9]
MWHernquist = MWProfile.HernquistMass(rarr,10,MWHaloMass[9])
'''
#plotting mass components and Hernquist for MW
plt.semilogy(rarr, MWHaloMass, color = 'red', label = 'Halo Mass')
plt.semilogy(rarr, MWDiskMass, color = 'blue', label = 'Disk Mass')
plt.semilogy(rarr, MWBuldgeMass, color = 'green', label = 'Buldge Mass')
plt.semilogy(rarr, MWTotMass, color = 'black', label = 'Total Mass')
plt.semilogy(rarr, MWHernquist, color = 'pink', label = 'Hernquist a=10')
legend = ax.legend(loc='lower right')
plt.title('MW Mass Profile')
plt.xlabel('Distance from Center of Mass (kpc)')
plt.ylabel('Mass Enclosed (Msun)')
#plt.show()

#Mass Components for M31
M31HaloMass = M31Profile.MassEnclosed(1,rarr)
M31DiskMass = M31Profile.MassEnclosed(2,rarr)
M31BuldgeMass = M31Profile.MassEnclosed(3,rarr)
M31TotMass = M31Profile.MassEnclosedTotal(rarr)
#Hernquist Mass Profile
M31Hernquist = M31Profile.HernquistMass(rarr,10,M31HaloMass[9])
#plotting mass components and Hernquist for MW
plt.semilogy(rarr, M31HaloMass, color = 'red', label = 'Halo Mass')
plt.semilogy(rarr, M31DiskMass, color = 'blue', label = 'Disk Mass')
plt.semilogy(rarr, M31BuldgeMass, color = 'green', label = 'Buldge Mass')
plt.semilogy(rarr, M31TotMass, color = 'black', label = 'Total Mass')
plt.semilogy(rarr, M31Hernquist, color = 'pink', label = 'Hernquist a=10')
legend = ax.legend(loc='lower right')
plt.title('M31 Mass Profile')
plt.xlabel('Distance from Center of Mass (kpc)')
plt.ylabel('Mass Enclosed (Msun)')
#plt.show()

#Mass Components for M33
M33HaloMass = M33Profile.MassEnclosed(1,rarr)
M33DiskMass = M33Profile.MassEnclosed(2,rarr)
M33TotMass = M33Profile.MassEnclosedTotal(rarr)
#Hernquist Mass Profile
M33Hernquist = M33Profile.HernquistMass(rarr,10.5,M33HaloMass[9])
#plotting mass components and Hernquist for M33
plt.semilogy(rarr, M33HaloMass, color = 'red', label = 'Halo Mass')
plt.semilogy(rarr, M33DiskMass, color = 'blue', label = 'Disk Mass')
plt.semilogy(rarr, M33TotMass, color = 'black', label = 'Total Mass')
plt.semilogy(rarr, M33Hernquist, color = 'pink', label = 'Hernquist a=10.5')
legend = ax.legend(loc='lower right')
plt.title('M33 Mass Profile')
plt.xlabel('Distance from Center of Mass (kpc)')
plt.ylabel('Mass Enclosed (Msun)')
#plt.show()
'''
#Velocity components of MW
vMWHalo = MWProfile.CircularVelocity(1,rarr)*(u.km/u.s)
print len(vMWHalo), type(vMWHalo)
vMWDisk = MWProfile.CircularVelocity(2,rarr)*(u.km/u.s)
vMWBuldge = MWProfile.CircularVelocity(3,rarr)*(u.km/u.s)
vMWTot = MWProfile.CircularVelocityTotal(rarr)*(u.km/u.s)
#Hernquist velocity
vMWHernquist = MWProfile.HernquistVCirc(rarr,1,MWHaloMass[9])
plt.plot(rarr, vMWHalo, color = 'red', label = 'Halo Velocity')
plt.plot(rarr, vMWDisk, color = 'blue', label = 'Disk Velocity')
plt.plot(rarr, vMWBuldge, color = 'green', label = 'Buldge Velocity')
plt.plot(rarr, vMWTot, color = 'black', label = 'Total Velocity')
#plt.plot(rarr, vMWHernquist, color = 'pink', label = 'Hernquist a=1')
legend = ax.legend(loc='lower right')
plt.title('MW Velocity Profile')
plt.xlabel('Distance from Center of Mass (kpc)')
plt.ylabel('Circular Velocity (km/s)')
plt.show()

#Velocity components of M31
vM31Halo = M31Profile.CircularVelocity(1,rarr)*(u.km/u.s)
#print len(vMWHalo), type(vMWHalo)
vM31Disk = M31Profile.CircularVelocity(2,rarr)*(u.km/u.s)
vM31Buldge = M31Profile.CircularVelocity(3,rarr)*(u.km/u.s)
vM31Tot = M31Profile.CircularVelocityTotal(rarr)*(u.km/u.s)
#Hernquist velocity
vM31Hernquist = M31Profile.HernquistVCirc(rarr,1,MWHaloMass[9])
plt.plot(rarr, vM31Halo, color = 'red', label = 'Halo Velocity')
plt.plot(rarr, vM31Disk, color = 'blue', label = 'Disk Velocity')
plt.plot(rarr, vM31Buldge, color = 'green', label = 'Buldge Velocity')
plt.plot(rarr, vM31Tot, color = 'black', label = 'Total Velocity')
#plt.plot(rarr, vM31Hernquist, color = 'pink', label = 'Hernquist a=1')
legend = ax.legend(loc='lower right')
plt.title('M31 Velocity Profile')
plt.xlabel('Distance from Center of Mass (kpc)')
plt.ylabel('Circular Velocity (km/s)')


           
        
            
            
            
   
            
            
            
        
