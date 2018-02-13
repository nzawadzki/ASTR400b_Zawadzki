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
        #use MassEnclosed to find the mass for each partice type at each radius
        mtot2 = 0
        count = 0
        n = np.arange(0,3)
        #looping over array element
        for element in range(len(r)):
            #find mass array of type n
            mtot1 = self.MassEnclosed(n,r)
            #add it to the total mass array
            mtotsum = np.sum(mtot1[element])
            #mtot2 is initialized as 0 before the loop in order for the addition to work for the first particle tyep
            mtotal.append(mtotsum)
            '''
            #caveat to account for M33 having no buldge
            if self.gname != 'M33' and n != 2:
                #reassign mtot2 for the next loop
                mtot2 = mtotal
            else:
                #do nothing
                a = 0'''
        
        mtotarray = mtotal
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

rarr = np.arange(0.1,10.1)

#Testing class
MWProfile = MassProfile('MW',0)
MWMass = MWProfile.MassEnclosed(2,rarr)
print MWMass

MWMassTot = MWProfile.MassEnclosedTotal(rarr)
print MWMassTot
    
           
        
            
            
            
   
            
            
            
        
