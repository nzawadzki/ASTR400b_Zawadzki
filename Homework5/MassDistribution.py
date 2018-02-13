#Determine the mass distribution of each galaxy at a certain snap shot and use this to determine each galaxy's rotation curve
import numpy as np
import astropy.units as u
#import necessary files
from ReadFile import Read
from CenterOfMass import CenterOfMass

#import relavent plotting modules
import matplotlib.pyplot as plt
#get_ipython().magic('matplotlib inline')

#constants
from astropy.constants import G

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
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        #Fg = G*Mm/r**2=mv**2/r
        #loop over radius array
        #define counter
        count = 0
        for radius in r:
            v = np.sqrt(G*self.MassEnclosed(ptype,r)[count]/radius)
            V.append(v)
            count = count + 1
            
        return round(V,2)*u.km/u.s
    
    #function to calculate the total circular velocity of (buldge+disk+halo) at each input radii
    def CircularVelocityTotal(self,r):
        return r
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
#radius array out to 30kpc
rarr = np.logspace(-1,1.477,10)
#Create class objects for each galaxy
MWProfile = MassProfile('MW',0)
M31Profile = MassProfile('M31',0)
M33Profile = MassProfile('M33',0)
#Mass Components of MW
MWHaloMass = MWProfile.MassEnclosed(1,rarr)
MWDiskMass = MWProfile.MassEnclosed(2,rarr)
MWBuldgeMass = MWProfile.MassEnclosed(3,rarr)
MWTotMass = MWProfile.MassEnclosedTotal(rarr)
#print len(MWHaloMass)

#Hernquist Mass Profile
MWHernquist = MWProfile.HernquistMass(rarr,10,MWHaloMass[9])
#plotting mass components and Hernquist
plt.semilogy(rarr, MWHaloMass, color = 'red', label = 'Halo Mass')
plt.semilogy(rarr, MWDiskMass, color = 'blue', label = 'Disk Mass')
plt.semilogy(rarr, MWBuldgeMass, color = 'green', label = 'Buldge Mass')
plt.semilogy(rarr, MWTotMass, color = 'black', label = 'Total Mass')
plt.semilogy(rarr, MWHernquist, color = 'pink', label = 'Hernquist a=10')
legend = ax.legend(loc='lower right')
plt.title('MW Mass Profile')
plt.xlabel('Distance from Center of Mass (kpc)')
plt.ylabel('Mass Enclosed (Msun)')
plt.show()

#MWMassTot = MWProfile.MassEnclosedTotal(rarr)
#print MWMassTot
    
           
        
            
            
            
   
            
            
            
        
